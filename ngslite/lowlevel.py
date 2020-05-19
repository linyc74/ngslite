import os
import subprocess
from typing import Optional


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


def _temp(prefix: str = 'temp', suffix: str = '') -> str:
    """
    Returns a temp file name that does not exist in the current folder, i.e. temp000.sam

    Args:
        prefix

        suffix:
            Usually the file extension
    """
    i = 0
    while True:
        name = f"{prefix}{i:03}{suffix}"
        if not os.path.exists(name):
            return name
        i += 1


def gzip(file: str, keep: bool = True):
    """
    Call the "gzip" command to zip or unzip files.

    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    is_gz = file.endswith('.gz')

    decompress = '--decompress ' if is_gz else ''

    if keep:
        stdout = '--stdout '
        if is_gz:
            output = f' > {file[:-3]}'
        else:
            output = f' > {file}.gz'
    else:
        stdout = ''
        output = ''

    call(f"gzip {decompress}{stdout}{file}{output}")


def printf(s: str):
    print(s, flush=True)
