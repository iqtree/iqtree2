from subprocess import PIPE, Popen, check_call

def exec_command(cmnd: str, stdout=PIPE, stderr=PIPE) -> str:
    """
    Executes a command and returns stdout if completes exit code 0.

    Parameters
    ----------
    cmnd : str
        The command to be executed.
    stdout, stderr : streams
        Default value (PIPE) intercepts process output, setting to None
        blocks this.

    Returns
    -------
    str
        The output of the command execution.
    """
    proc = Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        raise RuntimeError(f"FAILED: {cmnd}\n{msg}")
    return out.decode("utf8") if out is not None else None
