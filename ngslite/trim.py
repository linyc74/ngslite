from .lowlevel import __call
from functools import partial
printf = partial(print, flush=True)


def trim_galore(fq1, fq2, quality=20, gzip=True, length=20, log='trim_galore.log'):
    """
    Wrapper function of "trim_galore".

    Use the default settings of "trim_galore",
        e.g. Illumina P5 and P7 adapter sequences are used to trim reads.

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

        log: str, path-like
            Append the stderr of trim_galore to the <log> file
    """
    gzip = ['', '--gzip '][gzip]
    __call(f"trim_galore --paired --quality {quality} --phred33 --fastqc --illumina {gzip}--length {length} --max_n 0 --trim-n --retain_unpaired {fq1} {fq2} 2>> {log}")

