import os
import shutil
import unittest
from typing import Tuple


def get_dirs(py_path: str) -> Tuple[str, str, str]:
    indir = os.path.relpath(path=py_path[:-3], start=os.getcwd())
    basedir = os.path.dirname(indir)
    workdir = os.path.join(basedir, 'workdir')
    outdir = os.path.join(basedir, 'outdir')
    return indir, workdir, outdir


class TestCase(unittest.TestCase):

    def set_up(self, py_path: str):
        self.indir, self.workdir, self.outdir = get_dirs(py_path=py_path)
        for d in [self.workdir, self.outdir]:
            os.makedirs(d, exist_ok=True)

    def tear_down(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def assertFileEqual(self, file1: str, file2: str):
        with open(file1) as fh1:
            with open(file2) as fh2:
                self.assertEqual(fh1.read(), fh2.read())
