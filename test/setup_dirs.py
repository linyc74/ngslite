import os


def setup_dirs(fpath: str):
    """
    Args:
        fpath: The file path of the script given by __file__
    """
    cwd = os.getcwd()
    testdir = '.' if cwd.endswith('test') else './test'
    datadir = f'{testdir}/{os.path.basename(fpath)[:-3]}'
    workdir = f'{testdir}/workdir'
    outdir = f'{testdir}/outdir'

    os.makedirs(workdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    return testdir, datadir, workdir, outdir
