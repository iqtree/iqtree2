from subprocess import PIPE, Popen


def exec_command(cmnd: str, stdout=PIPE, stderr=PIPE) -> str:
    """
    Execute a command and capture the output.

    Parameters
    ----------
    cmnd : str
        The command to execute.
    stdout : int, optional
        The standard output for the command.
    stderr : int, optional
        The standard error for the command.

    Returns
    -------
    str
        The output of the command.
    """
    proc = Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        raise RuntimeError(f"FAILED: {cmnd}\n{msg}")
    return out.decode("utf8") if out is not None else None
