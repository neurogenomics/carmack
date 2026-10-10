import ctypes
import ctypes.util
import logging
import os
import subprocess
import sys
from functools import partial
from signal import SIGKILL

log = logging.getLogger(__name__)

LIBC = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)
PR_SET_PDEATHSIG = ctypes.c_int(1)  # <sys/prctl.h>


class LogSubprocess:
    """
    Helper class with functions for logging calls to subprocesses.
    """

    def __init__(self):
        """
        Initialise the LogSubprocess object
        """
        # Store the callback; prctl must never run in the calling process.
        self.pdeathsig = (
            self._child_preexec_set_pdeathsig if sys.platform.startswith("linux") else None
        )

    def _child_preexec_set_pdeathsig(self, parent_pid=None):
        """
        When used as the preexec_fn argument for subprocess.Popen etc,
        causes the subprocess to recieve SIGKILL if the parent process
        terminates.
        """
        if sys.platform.startswith("linux"):
            zero = ctypes.c_ulong(0)
            if LIBC.prctl(PR_SET_PDEATHSIG, ctypes.c_ulong(SIGKILL), zero, zero, zero) != 0:
                error = ctypes.get_errno()
                raise OSError(error, "Cannot set child parent-death signal")
            # The parent may have exited between fork and prctl. In that case
            # no later parent-death event will arrive, so close the race here.
            if parent_pid is not None and os.getppid() != parent_pid:
                os.kill(os.getpid(), SIGKILL)

    def _child_setup(self, parent_pid, caller_preexec):
        self.pdeathsig(parent_pid)
        if caller_preexec is not None:
            caller_preexec()

    def _with_child_setup(self, kwargs):
        if self.pdeathsig is not None:
            # Capture the actual launching PID, not the PID at wrapper creation:
            # a wrapper can be inherited by a multiprocessing worker.
            kwargs["preexec_fn"] = partial(
                self._child_setup, os.getpid(), kwargs.get("preexec_fn")
            )
        return kwargs

    def check_call(self, *args, **kwargs):
        return subprocess.check_call(*args, **self._with_child_setup(kwargs))

    def check_output(self, *args, **kwargs):
        return subprocess.check_output(*args, **self._with_child_setup(kwargs))

    def call(self, *args, **kwargs):
        return subprocess.call(*args, **self._with_child_setup(kwargs))

    def Popen(self, *args, **kwargs):
        return subprocess.Popen(*args, **self._with_child_setup(kwargs))
