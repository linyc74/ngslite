import subprocess


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def bowtie2_mapping(ref, fq1, sam, fq2=None):
    """
    Args:
        ref: str, path-like object
            The reference fasta

        fq1: str, path-like object
            The read-1 fastq

        fq2: str, path-like object
            The read-2 fastq. If none, use <fq1> for unpaired mapping.

        sam: str, path-like object
            The output SAM file
    """
    # Build the .bt2 index files
    cmd = 'bowtie2-build {} ref'.format(ref)
    __call(cmd)

    if fq2 is None:
        # Unpaired
        cmd = 'bowtie2 -x ref -U {} -S {}'.format(fq1, sam)
    else:
        # Paired-end
        cmd = 'bowtie2 -x ref -1 {} -2 {} -S {}'.format(fq1, fq2, sam)
    __call(cmd)

    # Remove the .bt2 index files
    cmd = 'rm *.bt2'
    __call(cmd)

