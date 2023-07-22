from subprocess import PIPE, Popen, check_call

def exec_command(cmnd: str, stdout=PIPE, stderr=PIPE) -> str:
    proc = Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        raise RuntimeError(f"FAILED: {cmnd}\n{msg}")
    return out.decode("utf8") if out is not None else None
