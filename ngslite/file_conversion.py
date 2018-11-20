import subprocess


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def sam_to_bam(file, keep=True):
    """
    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    # -S: input is a sam
    # -b: output is a bam
    # -h: include header section
    cmd = 'samtools view -S -b -h {} > {}.bam'.format(file, file[:-4])
    __call(cmd)
    if not keep:
        __call('rm {}'.format(file))


def bam_to_sam(file, keep=True):
    """
    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    # -h: include header section
    cmd = 'samtools view -h {} > {}.sam'.format(file, file[:-4])
    __call(cmd)
    if not keep:
        __call('rm {}'.format(file))


def fq_to_fa(file, keep=True):
    """
    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    if file.endswith('.fastq'):
        output = file[:-6] + '.fa'
    elif file.endswith('.fq'):
        output = file[:-3] + '.fa'
    cmd = 'seqtk seq -A {} > {}'.format(file, output)
    __call(cmd)
    if not keep:
        __call('rm {}'.format(file))

