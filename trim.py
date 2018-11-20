import subprocess


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def trim_galore(fq1, fq2, quality=20, gzip=True, length=20):
    """
    Args:
        fq1: str, path-like
            The read-1 fastq file

        fq2: str, path-like
            The read-2 fastq file

        quality: int
            phred33 score

        gzip: bool
            Compress the output fastq files or not

        length: int
            The minimal length (bp) of reads to be retained
    """
    gzip = ['', '--gzip '][gzip]
    __call(f"trim_galore --paired --quality {quality} --phred33 --illumina {gzip}--length {length} --max_n 0 --trim-n --retain_unpaired {fq1} {fq2}")

