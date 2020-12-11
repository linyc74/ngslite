from ..lowlevel import call


def trim_galore(
        fq1: str,
        fq2: str,
        quality: int = 20,
        gzip: bool = True,
        length: int = 20,
        log: str = 'trim_galore.log'):
    """
    Use the default settings of "trim_galore",
        e.g. Illumina P5 and P7 adapter sequences are used to trim reads

    Args:
        fq1:
            The read-1 fastq file

        fq2:
            The read-2 fastq file

        quality:
            phred33 score

        gzip:
            Compress the output fastq files or not

        length:
            The minimal length (bp) of reads to be retained

        log:
            Pipe stderr of trim_galore to the <log> file
    """
    gzip = ['', '--gzip '][gzip]
    call(f'trim_galore --paired --quality {quality} --phred33 --fastqc --illumina \
{gzip}--length {length} --max_n 0 --trim-n --retain_unpaired {fq1} {fq2} 2>> {log}')
