import os
import itertools
from typing import List, Optional
from .lowlevel import call, check_output


def get_files(
        source: str = '.',
        startswith: str = '',
        endswith: str = '',
        isfullpath: bool = False) -> List[str]:
    """
    Get all files that start with or end with some strings in the source folder

    Args:
        startswith

        endswith

        source: path-like
            The source directory

        isfullpath
            If True, return the full file paths in the list

    Returns:
        A list of file names, or file paths
        If no files found, return an empty list []
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


def get_dirs(
        startswith: str = '',
        endswith: str = '',
        source: str = '.',
        isfullpath: bool = False) -> List[str]:
    """
    Similar to 'get_files()' but finds directories
    """

    ret = []
    s, e = startswith, endswith
    for path, dirs, files in os.walk(source):
        if path == source:
            ret = [d for d in dirs if (d.startswith(s) and d.endswith(e))]

    if isfullpath:
        ret = [os.path.join(source, d) for d in ret]

    if ret:
        # Sort the file list so that the order will be
        #     consistent across different OS platforms
        ret.sort()
    return ret


def concat(files: List[str], output: str):
    """
    Args:
        files:
            Each str is a file path

        output: path-like
    """
    call(f"cat {' '.join(files)} > {output}")


def zip_broadcast(*lists) -> List[tuple]:
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


def tar(dir_: str,
        output: Optional[str] = None,
        verbose: bool = False) -> Optional[str]:
    """
    Args:
        dir_:
            Path to the input dir

        output:
            Path to the output .tar.gz file

        verbose:
            Verbose mode or not

    Returns:
        Path to the output .tar.gz file
    """

    # -c create archive
    # -z gzip
    # -v verbose
    # -f tar.gz filepath

    v = 'v' if verbose else ''

    if dir_.endswith('/'):
        dir_ = dir_[:-1]
    if output is None:
        output = f'{dir_}.tar.gz'
    f = f'-f "{output}"'

    current_dir = os.path.dirname(os.path.abspath(dir_))
    basename = os.path.basename(dir_)

    success = call(f'tar -cz{v} {f} -C "{current_dir}" "{basename}"')

    if success:
        return output
    else:
        return None


def untar(
        file: str,
        archive_path: Optional[str] = None,
        dstdir: Optional[str] = None,
        verbose: bool = False) -> Optional[str]:
    """
    Args:
        file:
            Path to the input .tar.gz file

        archive_path:
            Path to a specific archived file

        dstdir:
            Path to the destination directory in which files are extracted

        verbose:
            Verbose mode or not

    Returns:
        Path to the extracted directory, without the trailing '/'
    """

    # -x extract
    # -z gunzip
    # -v verbose
    # -f tar.gz filepath

    f = f'-f "{file}"'
    v = 'v' if verbose else ''
    a = '' if archive_path is None else f' "{archive_path}"'
    if dstdir is None:
        # dstdir is where the input file locates
        dstdir = os.path.dirname(file)

    success = call(f'tar -zx{v} {f} -C "{dstdir}"{a}')

    if not success:
        return None

    outpath = archive_path if archive_path else tar_list(file=file)[0]

    if outpath.endswith('/'):
        outpath = outpath[:-1]

    return os.path.join(dstdir, outpath)


def tar_list(file: str) -> List[str]:
    """
        Args:
            file:
                Path to the .tar.gz file

        Returns:
            List of archived directories and files
    """
    stdout = check_output(f'tar -tf "{file}"')
    return stdout.strip().split('\n')


def get_temp_path(
        prefix: str = 'temp',
        suffix: str = '') -> str:
    """
    Returns a temp file name that does not exist, i.e. temp000.txt

    Args:
        prefix:
            Can include the path

        suffix:
            Usually the file extension
    """
    i = 1
    while True:
        fpath = f'{prefix}{i:03}{suffix}'
        if not os.path.exists(fpath):
            return fpath
        i += 1


def gzip(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    gzip or gunzip files

    Args:
        file: path-like

        keep:
            Keep the input file or not

        output: path-like
    """
    is_gz = file.endswith('.gz')
    decompress = '--decompress' if is_gz else ''

    if output is None:
        if is_gz:
            output = f'{file[:-3]}'
        else:
            output = f'{file}.gz'

    args = [
        'gzip', decompress, '--stdout', f'"{file}"', ">", f'"{output}"'
    ]

    cmd = ' '.join(args)
    call(cmd=cmd)

    if not keep:
        os.remove(file)

    return output