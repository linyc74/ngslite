import subprocess


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def bedtools_multicov(bed, bams, output):
    """
    Args:
        bed: str, path-like
            The bed file, or any other interval file accepted by the bedtools

        bams: str, path-like; or list of str for multiple input files
            The bam file, or any other mapped read files accepted by the bedtools

        output: str, path-like
            The output file name
    """
    if isinstance(bams, list):
        bams = ' '.join(bams)

    cmd = 'bedtools multicov -bams {} -bed {} > {}'.format(bams, bed, output)
    __call(cmd)
