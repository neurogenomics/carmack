import csv
import logging
import os
from contextlib import ExitStack
from typing import Optional

import pysam

from carmack.bam_templates import PrimaryPairValidator, TemplateIdentity, template_identity
from carmack.utils import get_prefix, progress_bar

log = logging.getLogger(__name__)


class TagDedup:
    """
    Class that tags reads in BAM file with barcode and optionally removes
    duplicates.
    """

    def __init__(self, bam: str, bai: str, bc_valid_csv: str, umi_map: str | None = None) -> None:
        self.bam = bam
        self.bai = bai
        self.bc_valid_csv = bc_valid_csv
        self.umi_map = umi_map

        log.debug(
            f"TagDedup object created with BAM: {bam}, BAI: {bai},"
            f" valid barcodes CSV: {bc_valid_csv}, and UMI map: {umi_map}"
        )

    def _find_duplicate_templates(self, bam, bc_dict, umi_dict) -> set[TemplateIdentity]:
        """Choose a whole template per paired geometry; validate before writing.

        Read group scopes names, not molecule equivalence: duplicate comparisons
        still span lanes within the same cell. Unmapped/unpaired records are retained
        without coordinate-based collapse. Non-primary records inherit their primary
        template's decision and never choose a winner themselves.
        """
        with PrimaryPairValidator() as validator:
            seen_geometry: set[tuple] = set()
            duplicates: set[TemplateIdentity] = set()
            unpaired = unmapped = non_primary = missing_barcode = 0
            for read in bam.fetch(until_eof=True):
                if read.has_tag("CB"):
                    raise ValueError("Input BAM file reads already has 'CB' tags.")
                barcode = bc_dict.get(read.query_name)
                if barcode is None:
                    missing_barcode += 1
                    continue
                if read.is_secondary or read.is_supplementary:
                    non_primary += 1
                    continue
                if not read.is_paired:
                    unpaired += 1
                    continue
                identity = template_identity(read, barcode)
                identity = validator.observe(read, identity)
                if read.is_unmapped or read.mate_is_unmapped:
                    unmapped += 1
                    continue
                if not read.is_read1:
                    continue
                # Both ends are required. Sharing only one end is not evidence that two
                # complete fragments are duplicates; the decision is never made per mate.
                geometry = (
                    barcode,
                    read.reference_id,
                    read.reference_start,
                    read.is_reverse,
                    read.next_reference_id,
                    read.next_reference_start,
                    read.mate_is_reverse,
                    abs(read.template_length),
                )
                if umi_dict is not None:
                    geometry += (umi_dict.get(read.query_name, (None, None))[1],)
                if geometry in seen_geometry:
                    duplicates.add(identity)
                else:
                    seen_geometry.add(geometry)
            validated = validator.finish()
        log.info(
            f"Legacy dedup validated {validated} complete primary templates; "
            f"{len(duplicates)} duplicate templates. Record categories: "
            f"{unpaired} unpaired and {unmapped} unmapped-pair retained without collapse; "
            f"{non_primary} non-primary follow their template decision; "
            f"{missing_barcode} missing barcode skipped."
        )
        bam.reset()
        return duplicates

    def _tag(
        self,
        untagged_bam: pysam.AlignmentFile,
        tagged_bam: pysam.AlignmentFile,
        bc_dict: dict,
        umi_dict: dict[str, tuple[str, str]] | None = None,
        duplicate_templates: set[TemplateIdentity] | None = None,
    ) -> None:
        """
        Tag reads in untagged BAM file with barcode, duplicates and write to
        tagged BAM file.

        When ``umi_dict`` is provided, reads present in the map are additionally
        tagged with their raw (``UR``) and corrected (``UB``) UMI and the
        corrected UMI is used as an extra deduplication dimension. When
        ``umi_dict`` is ``None`` the behaviour is unchanged.
        """
        log.info(
            f"Starting to tag reads in BAM file. Total reads to process: "
            f"{untagged_bam.mapped + untagged_bam.unmapped}"
        )
        tag_count: int = 0

        duplicate_templates = duplicate_templates or set()
        total_reads = untagged_bam.mapped + untagged_bam.unmapped

        with progress_bar(unit="reads") as pbar:
            task = pbar.add_task("Tagging reads", total=total_reads)
            for read in untagged_bam.fetch(until_eof=True):
                read_name = read.query_name

                # Check if CB tag already exists
                if read.has_tag("CB"):
                    log.error(
                        f"Input BAM file read with CB tag detected (read name: {read_name}). Exiting."
                    )
                    raise ValueError("Input BAM file reads already has 'CB' tags.")

                # Log unpaired reads
                if not read.is_paired:
                    log.warning(f"Read {read_name} is not paired.")

                # Check if read has a barcode
                if read_name not in bc_dict:
                    log.warning(
                        f"Read {read_name} does not have a barcode in barcodes CSV"
                        " and will be skipped.",
                    )
                    continue

                # Get barcode and tag read
                barcode = bc_dict[read_name]
                read.set_tag("CB", barcode)
                # CR (raw cell barcode) currently mirrors CB: tag_dedup only
                # receives the corrected/matched barcode from bc_dict, so no
                # distinct raw cell barcode is available yet. Set CR to the same
                # value until a raw barcode is carried through the read->barcode
                # map (a future enhancement).
                read.set_tag("CR", barcode)
                log.debug(
                    f"Tagged read {read_name} (paired: {read.is_paired})"
                    f" with barcode {barcode}"
                )

                # Tag corrected UMI (UR/UB) when a UMI map is provided. Reads
                # absent from the map fall back to a None corrected UMI, so they
                # deduplicate exactly as in non-UMI mode.
                ub: str | None = None
                if umi_dict is not None and read_name in umi_dict:
                    ur, ub = umi_dict[read_name]
                    read.set_tag("UR", ur)
                    read.set_tag("UB", ub)

                read.set_tag(
                    "DU",
                    int(
                        read.is_paired and template_identity(read, barcode) in duplicate_templates
                    ),
                )

                # Write tagged read to tagged BAM file
                tagged_bam.write(read)
                tag_count += 1
                pbar.advance(task)

            log.info(f"Finished tagging reads. Total reads tagged:" f"{tag_count}")

    def _dedup(
        self,
        tagged_bam: pysam.AlignmentFile,
        dedup_bam: pysam.AlignmentFile,
        multiqc_log: Optional[csv.writer] = None,
    ):
        """
        Filter reads in tagged BAM file to remove duplicates.
        """
        log.info(
            f"Starting to deduplicate reads in tagged BAM file. Total "
            f"reads to process: {tagged_bam.mapped + tagged_bam.unmapped}"
        )

        unique_count: int = 0

        def get_dup_tag(read: pysam.AlignedSegment) -> bool:
            """
            Check and get duplicate tag from read.
            """
            if not read.has_tag("DU"):
                log.error("Read does not have a DU tag. Exiting.")
                raise ValueError("Read does not have a DU tag.")

            dup_tag = read.get_tag("DU")
            if not isinstance(dup_tag, int):
                log.error("Read has an invalid DU tag (non-int). Exiting.")
                raise ValueError("Read has an invalid DU tag (non-int).")

            if dup_tag not in [0, 1]:
                log.error("Read has an invalid DU tag (not 0 or 1). Exiting.")
                raise ValueError("Read has an invalid DU tag (not 0 or 1).")

            # Return True if duplicate, False if not
            return dup_tag == 1

        with progress_bar(unit="reads") as pbar:
            task = pbar.add_task(
                "Deduplicating reads", total=tagged_bam.mapped + tagged_bam.unmapped
            )
            for read in tagged_bam.fetch(until_eof=True):
                read_dup = get_dup_tag(read)

                if not read_dup:
                    dedup_bam.write(read)
                    unique_count += 1
                pbar.advance(task)

        # MultiQC log
        # Two reads for one read-pair
        if multiqc_log is not None:
            multiqc_log.writerow(["unique_reads", "duplicate_reads"])
            multiqc_log.writerow(
                [unique_count, tagged_bam.mapped + tagged_bam.unmapped - unique_count]
            )

    def tag_dedup_reads(self, dedup: bool, output_dir: str, prefix: Optional[str] = None) -> None:
        """
        Generate a tagged BAM file and TSV file with filtered reads coordinates
        with associated barcodes. If `dedup = TRUE`, also generate a
        deduplicated BAM file from the tagged BAM file and log duplicate counts.
        """
        # Set prefix, if non specified
        if prefix is None:
            prefix = get_prefix(self.bam)

        # Init
        BAM_TAGGED_PATH = os.path.join(output_dir, prefix + ".tagged.bam")
        BAM_DEDUP_PATH = os.path.join(output_dir, prefix + ".dedup.tagged.bam")
        MULTIQC_CSV_PATH = os.path.join(output_dir, prefix + ".dedup.stats_mqc.log")

        # Read barcodes from CSV
        with open(self.bc_valid_csv, "r") as valid_barcodes:
            csv_reader = csv.reader(valid_barcodes)
            bc_dict = {}
            for line in csv_reader:
                name, barcode = line[0].split(" ", 1)[0], line[1]
                if name in bc_dict and bc_dict[name] != barcode:
                    raise ValueError(f"Conflicting barcode assignments for read name {name!r}.")
                bc_dict[name] = barcode
        log.debug(f"Loaded {len(bc_dict)} valid barcodes.")

        # Read the corrected UMI map, if provided. Fixed upstream contract:
        # TAB-separated, no header, columns read_id, barcode, UR, UB.
        umi_dict: dict[str, tuple[str, str]] | None = None
        if self.umi_map is not None:
            with open(self.umi_map, "r") as umi_map_file:
                umi_reader = csv.reader(umi_map_file, delimiter="\t")
                umi_dict = {}
                for row in umi_reader:
                    name, barcode, ur, ub = row
                    if name in bc_dict and bc_dict[name] != barcode:
                        raise ValueError(
                            f"UMI-map barcode conflicts with barcode CSV for {name!r}."
                        )
                    if name in umi_dict and umi_dict[name] != (ur, ub):
                        raise ValueError(f"Conflicting UMI assignments for read name {name!r}.")
                    umi_dict[name] = (ur, ub)
            log.debug(f"Loaded {len(umi_dict)} UMI map entries.")

        # Tag and write to tagged BAM file
        with ExitStack() as stack:
            input_bam = stack.enter_context(
                pysam.AlignmentFile(self.bam, "rb", index_filename=self.bai)
            )
            duplicate_templates = self._find_duplicate_templates(input_bam, bc_dict, umi_dict)
            bam_tagged = stack.enter_context(
                pysam.AlignmentFile(BAM_TAGGED_PATH, "wb", header=input_bam.header)
            )

            # Tag
            self._tag(input_bam, bam_tagged, bc_dict, umi_dict, duplicate_templates)

        pysam.index(BAM_TAGGED_PATH)
        log.info(
            f"Tagged BAM file written to {BAM_TAGGED_PATH} with a corresponding BAI index file."
        )

        if not dedup:
            return

        # Deduplicate and write to deduplicated BAM file
        with ExitStack() as stack:
            bam_tagged = stack.enter_context(pysam.AlignmentFile(BAM_TAGGED_PATH, "rb"))
            bam_dedup = stack.enter_context(
                pysam.AlignmentFile(BAM_DEDUP_PATH, "wb", header=bam_tagged.header)
            )
            multiqc_csv = stack.enter_context(open(MULTIQC_CSV_PATH, mode="w", newline=""))
            multiqc_csv_writer = csv.writer(multiqc_csv)

            # Deduplicate
            self._dedup(bam_tagged, bam_dedup, multiqc_csv_writer)

        pysam.index(BAM_DEDUP_PATH)
        log.info(
            f"Deduplicated BAM file written to {BAM_DEDUP_PATH} with "
            "a corresponding BAI index file."
        )
        log.info(f"MultiQC log written to {MULTIQC_CSV_PATH}.")
