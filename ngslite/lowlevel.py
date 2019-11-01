from subprocess import check_call
from os.path import exists


def _call(cmd, print_cmd=True):
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


def _temp(prefix='temp', suffix=''):
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


def _gzip(file, keep=True):
    """
    Call the "gzip" command to zip or unzip files

    Args:
        file: str, path-like

        keep: bool, keep the input file or not
    """
    keep = ['', '-k '][keep]
    decomp = ['', '-d '][file.endswith('.gz')]
    _call(f"gzip {decomp}{keep}{file}")


def call(cmd, print_cmd=True):
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


def gzip(file, keep=True):
    """
    Call the "gzip" command to zip or unzip files.

    Args:
        file: str, path-like

        keep: bool, keep the input file or not
    """
    is_gz = file.endswith('.gz')

    decompress = '--decompress ' if is_gz else ''

    if keep:
        stdout = '--stdout '
        if is_gz: output = f' > {file[:-3]}'
        else: output = f' > {file}.gz'
    else:
        stdout = ''
        output = ''

    _call(f"gzip {decompress}{stdout}{file}{output}")


def printf(s):
    print(s, flush=True)
