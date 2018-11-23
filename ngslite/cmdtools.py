import os
import subprocess
import itertools


def __call(cmd, stderr=None):
    """
    This is a low-level method that takes a str command and execute it.
    When piping the stdout to a file, create a file object to be piped in
        by the subprocess.check_call() method.

    Args:
        cmd: str,
            The command to be executed

        stderr: str,
            The file name to write in standard error
    """
    # Open a stderr file if stderr is not None
    if stderr is not None:
        # stderr is usually written into a log file, so set append 'a' mode
        stderr = open(stderr, 'a')

    args = cmd.split(' ')

    # If the '>' is seen, then create a stdout file to pipe in
    if '>' in args:
        i = args.index('>')
        # Open a stdout file to which stdout is piped in
        stdout = open(args[i+1], 'w')
        args = args[:i]
    # If no stdout file, then set it to None
    # This will make the subprocess print out the stdout
    else:
        stdout = None

    print(f"CMD: {cmd}")
    try:
        subprocess.check_call(args, stdout=stdout, stderr=stderr)
    except Exception as inst:
        print(inst)

    if stdout is not None:
        stdout.close()
    if stderr is not None:
        stderr.close()


def __check_input_format(cmd, inputs):
    """
    This is a subprocess of the the method batch_call().
    Check the format of 'IN's in the <cmd> str and <inputs> (file names).
    """
    # num_IN: number of 'IN's in the cmd string
    num_IN = cmd.split(' ').count('IN')

    # If each input is a str, then there should only be one 'IN' in the cmd str
    if type(inputs[0]) == str:
        if num_IN != 1:
            print("[Error] There are {} 'IN's in the cmd but only 1 file for each input"
                  .format(num_IN)
                  )
            return False

    # If each input is a tuple, then len(each_input) should equal to number of 'IN's in the cmd str
    if type(inputs[0]) == tuple:
        if num_IN != len(inputs[0]):
            print("[Error] There are {} 'IN's in the cmd but {} files for each input".
                  format(num_IN, len(inputs[0])
                         )
                  )
            return False

    return True


def __check_output_format(cmd, outputs):
    """
    This is a subprocess of the method batch_call().
    Check the format of 'OUT' in the cmd string and outputs (file names)
    """
    # If cmd contains 'OUT' then the <output> arguments cannot be None
    if 'OUT' in cmd and outputs is None:
        print("[Error] <cmd> contains 'OUT' but <outputs> not specified.")
        return False
    return True


def call(cmd):
    """
    A simple wrapper of subprocess.check_call(<cmd>, shell=True).
    Stderr will be printed out in the console.

    Args:
        cmd: str
            The command to be executed
    """
    print(f"CMD: {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def batch_call(cmd, inputs, source='.', outputs=None, stderr=None):
    """
    A wrapper function of __call() to process an array of input files and output files.

    Input files could be in a specified <source> folder; output files is always in the current folder.

    An example of <inputs> and <outputs>:
        ['in_1.txt' , 'in_2.txt' , 'in_3.txt' ]
        ['out_1.txt', 'out_2.txt', 'out_3.txt']

    <inputs> and <outputs> must have the same length

    Args:
        cmd: str
            The command string. Use 'IN' and 'OUT' for the placeholder of
                input and output file names, respectively.
            When multiple input files are required for a single command,
                explicitly specify multiple 'IN' placeholders in the command str.

        inputs: list of str or tuple
            If each input is a single file, then each input element is a str.
            Sometimes each input is consist of multiple files,
                then use a tuple (of strings) to include multiple files.
            The input file list has to be of the same length of the output file list.

        source: str
            The path to the folder containing the input files. Default '.'

        outputs: list of str
            The output file list.

        stderr: str
            The file name to write in standard error
    """
    # Check the format of 'IN's in the <cmd> string and <inputs> (file names)
    if not __check_input_format(cmd, inputs):
        return

    # Check the format of 'OUT' in the <cmd> string and <outputs> (file names)
    if not __check_output_format(cmd, outputs):
        return

    # Check the lengths of <inputs> and <outputs>
    if outputs:
        if len(inputs) != len(outputs):
            print('[Error] Lengths of <inputs> and <outputs> are not the same.')
            return

    # For each input:
    #   (1) Replace 'IN's with input file paths
    #   (2) Replace the 'OUT' with output file path
    #   (3) Call the command
    for i, fnames in enumerate(inputs):
        args = cmd.split(' ')

        # If there's only one file for each input,
        #     convert the format to a tuple, as if there are multiple files for that input
        if type(fnames) == str:
            fnames = (fnames, )

        # Include the <source> file path to each file
        fnames = [os.path.join(source, f) for f in fnames]

        # Replace all 'IN's in the cmd with input files
        num_IN = cmd.split(' ').count('IN')
        for j in range(num_IN):
            args[args.index('IN')] = fnames[j]

        # For each command, there is going to be only one output file
        # So each output file is a single file, i.e. a str
        if 'OUT' in cmd:
            args[args.index('OUT')] = outputs[i]

        # Call the command
        __call(' '.join(args), stderr=stderr)


def iter_call(cmd, input_type, source='.', output_type=None, stderr=None):
    """
    A wrapper function of __call() to process all files in a folder with the same file extension <input_type>.
    The output files has the same filename as the input file, except with a new file extension <output_type>.

    Args:
        cmd: str,
            The command string. Use 'IN' and 'OUT' for the placeholder of
                input and output file names, respectively.

        input_type: str,
            Input file extension

        source: str,
            The path to the folder containing the input files. Default '.'

        output_type: str,
            Output file extension

        stderr: str,
            The file name to write in standard error
    """
    # Get a list of all input files with the extension <input_type>
    for path, dirs, files in os.walk(source):
        if path == source:
            inputs = [f for f in files if f.endswith(input_type)]
            inputs.sort()

    # For each input file
    for file in inputs:
        args = cmd.split(' ')

        # Replace 'IN' with the full path to the input file
        args[args.index('IN')] = os.path.join(source, file)
        if 'OUT' in args:
            # Replace 'OUT' with the output file
            # The output file extension is <output_type>
            args[args.index('OUT')] = file[:-len(input_type)] + output_type

        __call(' '.join(args), stderr=stderr)


def build_cmd(base_cmd, options, file_in=None, file_out=None):
    """
    Linux commands usually follow the format:
        <base_cmd> [options] [file_in] [> file_out]

        '<>' = required
        '[]' = optional

    This function builds a line of str of a command from the input arguments

    Args:
        base_cmd: str,
            The base command can contain blank space

        options: list of str or tuple, each tuple could be (str, str) or (str, int) or (str, float)
            The strings shouldn't contain any blank space.
            If a value is associated with an option, then use a tuple
            Here's an exmaple:
            ['-a', ('-b', 100)]

        file_in:
            str, for single file. 'IN' can be used for the placeholder of input files
            or None if there's no input

        file_out:
            str, for single. 'OUT' can be used for the placeholder of output files
            or None if there's no output

    Returns:
        str of the command
    """
    cmd = base_cmd
    for o in options:
        if type(o) == str:
            cmd = cmd + ' ' + o
        elif type(o) in (tuple, list):
            # Sometimes each option 'o' is a tuple with (--option, value), with the value being an integer or float
            # So convert the value to str by using map(str, o)
            cmd = cmd + ' ' + ' '.join(map(str, o))
    if file_in is not None:
        cmd = cmd + ' ' + file_in
    if file_out is not None:
        cmd = cmd + ' > ' + file_out
    return cmd

