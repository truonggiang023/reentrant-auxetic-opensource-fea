"""Interrupt real nested process groups without launching an FEA solver."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

from calibration_study import ROOT, prepare


def running(pid):
    status = Path(f'/proc/{pid}/stat')
    if not status.exists():
        return False
    try:
        return status.read_text().split()[2] != 'Z'
    except FileNotFoundError:
        return False


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux process-group cleanup contract')
def test_interrupting_driver_stops_separately_detached_native_worker(tmp_path):
    worker_pid, cli_pid = tmp_path/'worker.pid', tmp_path/'cli.pid'
    worker = ("import os,signal,time; from pathlib import Path; "
              "signal.signal(signal.SIGINT,signal.SIG_IGN); signal.signal(signal.SIGTERM,signal.SIG_IGN); "
              f"Path({str(worker_pid)!r}).write_text(str(os.getpid())); time.sleep(60)")
    child = ("import os,sys; from pathlib import Path; from run_radioss import run_logged; "
             f"Path({str(cli_pid)!r}).write_text(str(os.getpid())); "
             f"run_logged([sys.executable,'-c',{worker!r}], {str(ROOT)!r}, {str(tmp_path/'native.log')!r},60)")
    driver = ("import os,sys; from run_calibration import run_logged\n"
              "try:\n"
              f" run_logged([sys.executable,'-c',{child!r}],{str(tmp_path/'driver.log')!r},os.environ.copy(),60)\n"
              "except KeyboardInterrupt:\n sys.exit(130)\n")
    process = subprocess.Popen([sys.executable, '-c', driver], cwd=ROOT,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pids = []
    try:
        deadline = time.monotonic()+10
        while not worker_pid.exists() and time.monotonic() < deadline:
            if process.poll() is not None:
                stdout, stderr = process.communicate()
                pytest.fail(f'Driver exited early: {stdout!r} {stderr!r}')
            time.sleep(.02)
        assert worker_pid.exists(), 'Detached worker did not start'
        pids = [int(cli_pid.read_text()), int(worker_pid.read_text())]
        process.send_signal(signal.SIGINT)
        assert process.wait(timeout=18) == 130
        assert not any(running(pid) for pid in pids)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=3)
        for file in (cli_pid, worker_pid):
            if file.exists():
                try:
                    os.killpg(int(file.read_text()), signal.SIGKILL)
                except ProcessLookupError:
                    pass


def test_interrupted_execution_preserves_partial_output_journals_and_exits_130(tmp_path, monkeypatch):
    import run_calibration
    manifest = prepare(ROOT/'configs/dong2019_brick_damage.json', tmp_path/'study')
    row = next(row for row in manifest['runs']
               if row['candidate'] == 'locus_0_1' and row['role'] == 'training')
    def interrupt(*args):
        deck = Path(row['deck'])
        deck.parent.mkdir(parents=True, exist_ok=True)
        deck.write_text('Partial solver evidence must survive interruption\n')
        raise KeyboardInterrupt
    monkeypatch.setattr(run_calibration, 'run_logged', interrupt)
    monkeypatch.setattr(sys, 'argv', ['run_calibration.py', str(tmp_path/'study/manifest.json'),
                                    '--candidate', 'locus_0_1', '--role', 'training'])
    assert run_calibration.main() == 130
    journal = json.loads((tmp_path/'study/execution-journal.json').read_text())
    assert journal['events'][-1]['status'] == 'interrupted_preserved'
    assert Path(row['deck']).read_text() == 'Partial solver evidence must survive interruption\n'
