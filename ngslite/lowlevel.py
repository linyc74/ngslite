from subprocess import check_call
from os.path import exists


def __call(cmd, print_cmd=True):
    """
    Args:
        cmd: str
        print_cmd: bool
    """
    if print_cmd:
        print(f"CMD: {cmd}", flush=True)
    try:
        check_call(cmd, shell=True)
    except Exception as inst:
        print(inst, flush=True)


def __temp(prefix='temp', suffix=''):
    """
    Returns a temp file name that does not exists in the current folder, i.e. temp000.sam

    Args:
        prefix: str

        suffix: str
            Usually for the file extension

    Returns: str
    """
    i = 0
    while True:
        name = f"{prefix}{i:03}{suffix}"
        if not exists(name):
            return name
        i += 1
