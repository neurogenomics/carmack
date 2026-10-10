"""Child-only setup, caller callback composition, and Linux lifecycle contracts."""

import ctypes
import os
import sys
from signal import SIGKILL
from unittest.mock import MagicMock

import pytest

from carmack.io import log_subprocess
from tests.linux_parent_death_probe import probe


@pytest.fixture
def linux(monkeypatch):
    libc = MagicMock()
    libc.prctl.return_value = 0
    monkeypatch.setattr(log_subprocess, "LIBC", libc)
    monkeypatch.setattr(log_subprocess.sys, "platform", "linux")
    return libc


@pytest.mark.parametrize("method", ["Popen", "call", "check_call", "check_output"])
def test_callbacks_are_deferred_and_composed(monkeypatch, linux, method):
    invoked = []
    wrapper = log_subprocess.LogSubprocess()
    monkeypatch.setattr(log_subprocess.os, "getpid", lambda: 123)
    monkeypatch.setattr(log_subprocess.os, "getppid", lambda: 123)
    launcher = MagicMock()
    monkeypatch.setattr(log_subprocess.subprocess, method, launcher)
    getattr(wrapper, method)(["cmd"], preexec_fn=lambda: invoked.append("caller"))
    linux.prctl.assert_not_called()
    assert invoked == []
    launcher.call_args.kwargs["preexec_fn"]()
    linux.prctl.assert_called_once()
    assert invoked == ["caller"]


def test_parent_exits_before_prctl_race_kills_child(monkeypatch, linux):
    monkeypatch.setattr(log_subprocess.os, "getppid", lambda: 1)
    monkeypatch.setattr(log_subprocess.os, "getpid", lambda: 456)
    kill = MagicMock()
    monkeypatch.setattr(log_subprocess.os, "kill", kill)
    log_subprocess.LogSubprocess()._child_preexec_set_pdeathsig(parent_pid=123)
    kill.assert_called_once_with(456, SIGKILL)


def test_prctl_failure_is_not_ignored(linux):
    linux.prctl.return_value = -1
    ctypes.set_errno(22)
    with pytest.raises(OSError, match="parent-death signal"):
        log_subprocess.LogSubprocess()._child_preexec_set_pdeathsig(parent_pid=os.getpid())


def test_non_linux_preserves_caller_preexec(monkeypatch):
    monkeypatch.setattr(log_subprocess.sys, "platform", "darwin")
    launcher, caller = MagicMock(), MagicMock()
    monkeypatch.setattr(log_subprocess.subprocess, "Popen", launcher)
    log_subprocess.LogSubprocess().Popen(["cmd"], preexec_fn=caller)
    assert launcher.call_args.kwargs["preexec_fn"] is caller
    caller.assert_not_called()


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Linux prctl lifecycle")
def test_actual_linux_parent_death_and_session():
    assert probe()["parent_signal_unchanged"]
