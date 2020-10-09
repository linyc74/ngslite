import os
import unittest
from typing import Tuple


def setup_dirs(fpath: str) -> Tuple[str, str, str]:
    """
    Args:
        fpath: The file path of the script given by __file__
    """
    indir = fpath[:-3]
    workdir = f'{os.path.dirname(indir)}/workdir'
    outdir = f'{os.path.dirname(indir)}/outdir'

    os.makedirs(workdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    return indir, workdir, outdir


class TestCase(unittest.TestCase):

    def assertFileEqual(self, file1: str, file2: str):
        with open(file1) as fh1:
            with open(file2) as fh2:
                self.assertEqual(fh1.read(), fh2.read())
