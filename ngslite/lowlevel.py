import subprocess
from typing import Optional, Any


def call(cmd: str, print_cmd: bool = True):

    if print_cmd:
        print(f"CMD: {cmd}", flush=True)

    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst, flush=True)


def check_output(cmd: str) -> Optional[str]:
    """
    Wrapper function of subprocess.check_output

    This is meant to act like a function that returns the stdout of a command,
        but not a method that does something, so do not log the cmd

    Returns:
        stdout
        None if call error
    """
    try:
        stdout = subprocess.check_output(cmd, shell=True)
        return str(stdout, 'utf-8')

    except Exception as inst:
        print(inst, flush=True)
        return None


def printf(s: Any):
    print(s, flush=True)
