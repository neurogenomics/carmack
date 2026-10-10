"""Tiny Linux-only integration probe (standard library + Carmack, no test dependency).

Run from the repository root: python -m tests.linux_parent_death_probe
Children are short-lived, uniquely identified, and killed on a failed assertion.
"""

import ctypes
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from carmack.io.log_subprocess import LogSubprocess
from carmack.io.subprocess_stream import SubprocessStream

QUERY = """
import ctypes,json,os
libc=ctypes.CDLL(None)
value=ctypes.c_int()
assert libc.prctl(2,ctypes.byref(value),0,0,0)==0
print(json.dumps(dict(signal=value.value,pid=os.getpid(),sid=os.getsid(0))),flush=True)
"""


def probe():
    if not sys.platform.startswith("linux"):
        raise RuntimeError("This lifecycle probe requires Linux")
    libc = ctypes.CDLL(None)
    before = ctypes.c_int()
    assert libc.prctl(2, ctypes.byref(before), 0, 0, 0) == 0
    wrapper = LogSubprocess()
    child = json.loads(wrapper.check_output([sys.executable, "-c", QUERY]))
    assert child["signal"] == signal.SIGKILL
    # A custom preexec callback must compose with, not replace, parent-death setup.
    child_with_callback = json.loads(
        wrapper.check_output([sys.executable, "-c", QUERY], preexec_fn=os.setsid)
    )
    assert child_with_callback["signal"] == signal.SIGKILL
    assert child_with_callback["pid"] == child_with_callback["sid"]
    with SubprocessStream([sys.executable, "-c", QUERY]) as stream:
        streamed = json.loads(b"".join(stream))
    assert streamed["signal"] == signal.SIGKILL
    assert streamed["sid"] == streamed["pid"]
    after = ctypes.c_int()
    assert libc.prctl(2, ctypes.byref(after), 0, 0, 0) == 0
    assert before.value == after.value

    results = []
    with tempfile.TemporaryDirectory(prefix="carmack-pdeath-") as directory:
        for stream_mode in (False, True):
            marker = Path(directory) / ("stream" if stream_mode else "direct")
            child_code = (
                "import os,time; from pathlib import Path; "
                f"Path({str(marker)!r}).write_text(str(os.getpid())); time.sleep(30)"
            )
            parent_code = "\n".join(
                [
                    "import os,sys,time",
                    "from pathlib import Path",
                    "from carmack.io.log_subprocess import LogSubprocess",
                    "from carmack.io.subprocess_stream import SubprocessStream",
                    f"cmd=[sys.executable,'-c',{child_code!r}]",
                    (
                        "child=SubprocessStream(cmd)"
                        if stream_mode
                        else "child=LogSubprocess().Popen(cmd)"
                    ),
                    f"marker=Path({str(marker)!r})",
                    "deadline=time.monotonic()+5",
                    "while not marker.exists() and time.monotonic()<deadline: time.sleep(0.01)",
                    "os._exit(0 if marker.exists() else 1)",
                ]
            )
            parent = subprocess.Popen([sys.executable, "-c", parent_code])
            pid = None
            exited = False
            try:
                assert parent.wait(timeout=10) == 0
                pid = int(marker.read_text())
                deadline = time.monotonic() + 5
                state = None
                while time.monotonic() < deadline:
                    try:
                        state = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[0]
                    except FileNotFoundError:
                        state = "gone"
                    # A zombie has exited; reaping belongs to the host's init/subreaper.
                    if state in ("Z", "gone"):
                        break
                    time.sleep(0.01)
                assert state in ("Z", "gone"), (pid, state)
                exited = True
                results.append(
                    {"mode": "stream" if stream_mode else "direct", "child_state": state}
                )
            finally:
                if parent.poll() is None:
                    parent.kill()
                    parent.wait()
                if pid is not None and not exited:
                    try:
                        os.kill(pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
    result = {
        "parent_signal_unchanged": True,
        "child_signal": child["signal"],
        "stream_has_own_session": True,
        "parent_exit_checks": results,
    }
    print(json.dumps(result, sort_keys=True))
    return result


if __name__ == "__main__":
    probe()
