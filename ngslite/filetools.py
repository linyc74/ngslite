from .lowlevel import __call
import os
import itertools
from functools import partial
printf = partial(print, flush=True)


def zip_broadcast(*lists):
    """
    This is a broadcasting version of the built-in zip() function
    """
    list_of_lists = []

    # Get the longest length among the input lists
    longest = max(map(len, lists))

    for L in lists:
        if len(L) < longest:
            # If L is not the longest list, use itertools to cycle (i.e. broadcast) the list
            list_of_lists.append(itertools.cycle(L))
        else:  # len(L) == longest
            # If L is the longest list, than just use L to zip
            # len(L) defines the total length of the zipped result
            list_of_lists.append(list(L))

    return list(zip(*list_of_lists))


def get_files(startswith='', endswith='', source='.', isfullpath=False):
    """
    Get all files that start with or end with some strings in the source folder.

    Args:
        startswith: str

        endswith: str

        source: str, path-like
            The source directory.

        isfullpath: bool,
            If True, return the full file paths in the list.

    Returns:
        A list of file names, or file paths.
        If no files found, return an empty list [].
    """
    ret = []
    s, e = startswith, endswith
    for path, dirs, files in os.walk(source):
        if path == source:
            ret = [f for f in files if (f.startswith(s) and f.endswith(e))]

    if isfullpath:
        ret = [os.path.join(source, f) for f in ret]

    if ret:
        # Sort the file list so that the order will be
        #     consistent across different OS platforms
        ret.sort()
    return ret


def change_extension(files, old, new):
    """
    Change suffix of a list of file names.

    If a file does not have the specified old suffix, then it is not included
    in the returned list.

    Args
        files: list of str
        old: str, the old suffix
        new: str, the new suffix

    Returns:
        A list of file names with new suffices
    """
    ret = []
    for f in files:
        if f.endswith(old):
            ret.append(f[:-len(old)] + new)
    return ret


def change_prefix(files, old, new):
    """
    Change prefix of a list of file names.

    If a file does not have the specified old suffix, then it is not included
    in the returned list.

    Args
        files: list of str
        old: str, the old prefix
        new: str, the new prefix

    Returns:
        A list of file names with new suffices
    """
    ret = []
    for f in files:
        if f.startswith(old):
            ret.append(new + f[len(old):])
    return ret


def concat(files, output):
    """
    Args:
        files: list (or tuple) of str
            Each str is a file path

        output: str, path-like
    """
    __call(f"cat {' '.join(files)} > {output}")


def gzip(file, keep=True):
    """
    Call the "gzip" command to zip or unzip files.

    Args:
        file: str, path-like

        keep: bool, keep the input file or not
    """
    keep = ['', '-k '][keep]
    decomp = ['', '-d '][file.endswith('.gz')]
    __call(f"gzip {decomp}{keep}{file}")


def call(cmd):
    """
    A simple wrapper of subprocess.check_call(<cmd>, shell=True)

    Args:
        cmd: str
            The command to be executed
    """
    __call(cmd)
