import subprocess


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


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


def get_files(endswith='', startswith='', source='.'):
    """
    Get all files endswith a certain extension in the source folder

    Args:
        endswith: str
        source: str

    Returns:
        A list of file names. Just file names, not the full path.
        If no files found, return an empty list [].
    """
    ret = []
    s, e = startswith, endswith
    for path, dirs, files in os.walk(source):
        if path == source:
            ret = [f for f in files if (f.startswith(s) and f.endswith(e))]
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
    files = ' '.join(files)
    cmd = 'cat {} > {}'.format(files, output)
    __call(cmd)

