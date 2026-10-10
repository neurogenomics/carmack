import csv
import logging
import os
from typing import Dict, Generator, Optional, Tuple

import matplotlib.pyplot as plt
import pysam
from kneed import KneeLocator
from scipy.io import mmwrite

from carmack import __version__
from carmack.cell_caller.peak_barcode_matrix import PeakBarcodeMatrix
from carmack.io.bed_file import BedFile
from carmack.utils import get_prefix

log = logging.getLogger(__name__)


class CellCaller:
    """
    Class that calculates peak-barcode matrix from tagged BAM alignment file and peak file. Cells
    can then be filtered out based on the number of peaks they have dynamically.
    """

    def __init__(self, bed: str, bam: str, bai: Optional[str] = None) -> None:
        self.bam_path = bam
        self.bai_path = bai
        self.bed_path = bed

        self.matrix = None  # Later initialised as PeakBarcodeMatrix
        self.peaks_by_name = {}

        log.debug(f"CellCaller object created with BAM: {bam}, BAI: {bai}, and BED: {bed}")

    def check_matrix_exists(self, force: bool = False) -> None:
        """
        Check if matrix is already initialised. Raise error, if not.
        """
        if self.matrix is not None:
            if force:
                log.warning("Matrix already initialised. Overwriting...")
            else:
                raise AttributeError("Matrix already initialised. Use force=True to overwrite.")

    def get_overlaps(
        self, peak: Dict, bam: pysam.AlignmentFile, min_overlap: int
    ) -> Generator[str, None, None]:
        """
        Return all barcodes from reads which overlap with the input read. Only consider peaks that
        overlap with at least `min_overlap` bases.
        """
        chrom = peak["chrom"]
        start = peak["start"]
        end = peak["end"]

        for overlap in bam.fetch(chrom, start, end):
            # Count aligned query bases, excluding deletions and skipped reference regions.
            overlap_bases = overlap.get_overlap(start, end) or 0
            if overlap_bases >= min_overlap:
                yield overlap.get_tag("CB")

    def compute_matrix(self, min_overlap: int = 1, force: bool = False) -> None:
        """
        Find all overlapping peaks (with min_overlap bases) for each parcode and fill the
        peak-barcode matrix.

        The resultant matrix is stored in the `matrix` attribute.
        """
        log.debug("Computing peak-barcode matrix...")

        self.check_matrix_exists(force)

        if min_overlap < 1:
            raise ValueError("Minimum overlap should be at least 1 base.")

        # Initialise readers and matrix
        bed = BedFile(self.bed_path)
        peaks = []
        seen_entries = set()
        for entry in bed.open_read_iterator():
            key = (entry["chrom"], entry["start"], entry["end"])
            if key not in seen_entries:
                seen_entries.add(key)
                peaks.append(entry)
        peaks = tuple(peaks)
        peak_names = [peak["name"] for peak in peaks]
        if len(set(peak_names)) != len(peak_names):
            raise ValueError(
                "Duplicate peak names identify different intervals. Assign a unique BED "
                "name to each interval before calling cells."
            )
        with pysam.AlignmentFile(self.bam_path, "rb", index_filename=self.bai_path) as bam:
            barcodes = {read.get_tag("CB") for read in bam}
        self.matrix = PeakBarcodeMatrix(peak_names=peak_names, barcodes=barcodes, force=force)
        # Preserve exactly the intervals counted, even if the source BED changes before export.
        self.peaks_by_name = {peak["name"]: peak for peak in peaks}

        # Compute matrix
        with pysam.AlignmentFile(self.bam_path, "rb", index_filename=self.bai_path) as bam:
            for peak in peaks:
                peak_name = peak["name"]

                overlap_barcodes = self.get_overlaps(peak, bam, min_overlap)

                for barcode in overlap_barcodes:
                    self.matrix.increment_index(peak_name, barcode)

        log.debug(f"Peak-barcode matrix computed with {min_overlap}bp of minimum overlap.")

    def find_knee(self) -> Tuple[int, int]:
        """
        Find the knee point in the barcode rank plot. Uses the KneeLocator class from the kneed
        package, with the curve set to "concave" and direction set to "decreasing".
        """
        if self.matrix is None:
            raise AttributeError(
                "Matrix not computed. Run `compute_matrix` before computing knee."
            )

        barcode_sums = sorted(self.matrix.sum_barcodes(), reverse=True)
        if len(barcode_sums) < 2 or barcode_sums[0] == barcode_sums[-1]:
            raise ValueError(
                "No knee can be determined from an empty, singleton or flat barcode-rank "
                "curve. Inspect the counts and supply force_n explicitly (0 selects all)."
            )
        ranks = list(range(1, len(barcode_sums) + 1))

        knee_locator = KneeLocator(
            ranks, barcode_sums, curve="concave", direction="decreasing", online=True
        )
        knee_point_x = knee_locator.knee
        knee_point_y = knee_locator.knee_y
        if knee_point_x is None or knee_point_y is None:
            raise ValueError(
                "No knee was detected in the barcode-rank curve. Inspect the counts and "
                "supply force_n explicitly (0 selects all)."
            )

        log.debug(f"Knee point x (n_cells): {knee_point_x}")
        log.debug(f"Knee point y (minimum alignment overlaps per barcode): {knee_point_y}")

        return knee_point_x, knee_point_y

    def check_n(self, force_n: int) -> int:
        """
        Check if the number of cells to be selected is valid.
        """
        if not isinstance(force_n, int):
            raise TypeError("Number of cells to be selected should be an integer.")
        if force_n > len(self.matrix.barcodes):
            raise ValueError(
                "Number of cells to be selected is greater than the total number of cells. "
                'Use "0" to select all cells.'
            )
        if force_n < 0:
            raise ValueError(
                'Number of cells to be selected cannot be negative. Use "0" to select all cells.'
            )
        log.debug(f"Forcing selection of top {force_n} cells.")

        if force_n == 0:
            log.info("Forcing selection of all cells.")
            return len(self.matrix.barcodes)

        return force_n

    def make_plot(self, force_n: int = None) -> plt.Figure:
        """
        Create a barcode rank plot to visualise the number of overlapping peaks per barcode.
        Also highlights the cells to be selected based on the knee point or forced selection.
        """
        if self.matrix is None:
            raise AttributeError("Matrix not computed. Run `compute_matrix` before plotting.")

        log.debug("Creating barcode rank plot...")

        barcode_sums = sorted(self.matrix.sum_barcodes(), reverse=True)
        if not barcode_sums:
            raise ValueError("Cannot plot an empty barcode-rank curve.")
        ranks = list(range(1, len(barcode_sums) + 1))

        if force_n is not None:
            # Force n-cells
            log.debug(f"Forcing selection of top n-cells. ({force_n})")
            thresh_x = self.check_n(force_n)
            thresh_y = barcode_sums[thresh_x - 1]
            mode = "force_n_cells"
        else:
            # Use knee point
            log.debug("Finding knee point in barcode rank plot...")
            thresh_x, thresh_y = self.find_knee()
            mode = "knee_point"

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(ranks[:thresh_x], barcode_sums[:thresh_x], color="orange", label="Selected cells")
        ax.plot(
            ranks[thresh_x - 1 :],
            barcode_sums[thresh_x - 1 :],
            color="blue",
            label="Excluded cells",
        )

        # Threshold labels
        if force_n != 0:
            ax.axvline(thresh_x, color="gray", alpha=0.5)
            ax.text(
                x=thresh_x + 1,
                y=0.96,
                s=f"{thresh_x} cells",
                transform=ax.get_xaxis_transform(),
                ha="left",
                va="top",
            )
            ax.axhline(thresh_y, color="gray", alpha=0.5)
            ax.text(
                x=0.96,
                y=thresh_y + 1,
                s=f"{thresh_y} alignment overlaps",
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="bottom",
            )

        ax.text(
            y=-0.09,
            x=max(ranks),
            s=f"Mode: {mode}",
            transform=ax.get_xaxis_transform(),
            ha="right",
            va="bottom",
        )

        ax.set_xscale("log")
        ax.set_yscale("log")

        # Customise more
        ax.legend(loc="upper right")
        ax.grid(True, which="both", ls="-", alpha=0.2)
        ax.set_xlabel("Barcodes")
        ax.set_ylabel("Alignment overlaps with peaks")
        ax.set_title("Barcode Rank Plot")
        ax.set_xlim(left=1)
        fig.tight_layout()

        log.info("Barcode rank plot created.")

        return fig

    def export(self, output_dir: str, prefix: str, force_n: int = None) -> None:
        """
        Export to standard single-cell formats. If `force_n` is supplied, only the top n cells
        are returned. Otherwise, the knee point is used to determine the number of cells to be
        selected.

        Output files: barcodes.tsv, peaks.bed and matrix.mtx (with prefix, if supplied)

        Additional file description: Cell Ranger ATAC count: Feature-Barcode Matrices
        https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices
        """
        log.info("Starting export...")
        if self.matrix is None:
            raise AttributeError("Matrix not computed. Run `compute_matrix` before exporting.")

        # Prefix/output file path setup
        if prefix is None:
            prefix = get_prefix(self.bam_path)

        BARCODES_TSV_PATH = os.path.join(output_dir, f"{prefix}_barcodes.tsv")
        PEAKS_BED_PATH = os.path.join(output_dir, f"{prefix}_peaks.bed")
        MATRIX_MTX_PATH = os.path.join(output_dir, f"{prefix}_matrix.mtx")

        # Cell selection
        if force_n is None:
            n, _ = self.find_knee()
        elif force_n == 0:
            n = len(self.matrix.barcodes)
        else:
            n = self.check_n(force_n)

        log.info(f"Exporting top {n} cells.")

        # Filter barcodes
        barcode_with_counts = zip(self.matrix.barcodes, self.matrix.sum_barcodes())
        selected_barcodes = sorted(barcode_with_counts, key=lambda x: x[1], reverse=True)[:n]
        min_fragments = min((count for _, count in selected_barcodes), default=0)  # Logging
        selected_barcodes = [barcode for barcode, _ in selected_barcodes]

        log.debug(f"Minimum alignment overlaps per barcode in selected cells: {min_fragments}")

        # Filter peaks
        # Only return peaks that have at least one overlap with the selected barcodes, reduces size
        selected_peaks = [
            peak
            for peak in self.matrix.peaks
            if any(self.matrix.get_value(peak, barcode) for barcode in selected_barcodes)
        ]

        # Filter matrix
        # Only return the selected barcodes and peaks
        selected_barcodes_indices = [
            self.matrix.barcodes.index(barcode) for barcode in selected_barcodes
        ]
        selected_peaks_indices = [self.matrix.peaks.index(peak) for peak in selected_peaks]
        selected_matrix = self.matrix.filter_matrix(
            selected_peaks_indices, selected_barcodes_indices
        )

        log.debug(
            f"Shape of the matrix to be exported (peaks x barcodes): {selected_matrix.shape}"
        )

        # Export barcodes
        with open(BARCODES_TSV_PATH, "w") as f:
            writer = csv.writer(f, delimiter="\t")
            for barcode in selected_barcodes:
                writer.writerow([barcode])

        log.info(f"Exported barcodes to {BARCODES_TSV_PATH}")

        # Export peaks
        with open(PEAKS_BED_PATH, "w") as f:
            writer = csv.writer(f, delimiter="\t")
            # Matrix rows follow self.matrix.peaks, which is sorted by name, not BED order.
            for peak_name in selected_peaks:
                peak = self.peaks_by_name[peak_name]
                writer.writerow([peak["chrom"], peak["start"], peak["end"]])

        log.info(f"Exported peaks to {PEAKS_BED_PATH}")

        # Export matrix
        mmwrite(
            MATRIX_MTX_PATH,
            selected_matrix,
            comment=f"Peak-barcode matrix prepared with Carmack (v{__version__})",
        )

        log.info(f"Exported matrix to {MATRIX_MTX_PATH} (shape={selected_matrix.shape})")
        log.info("Export complete.")
