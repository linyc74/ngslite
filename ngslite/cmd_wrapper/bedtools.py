from typing import Union, List
from ..lowlevel import call
from ..filetools import get_temp_path


def bedtools_multicov(
        bed: str,
        bams: Union[str, List[str]],
        output: str):
    """
    Wrapper function of the command "bedtools multicov -bams <bams> -bed <bed> > <output>"
    Adds a header line to the <output> file

    Args:
        bed:
            The bed file path, or any other interval file accepted by the bedtools

        bams:
            One path or a list of paths
            The bam file, or any other mapped read files accepted by the bedtools

        output:
            The output tab-separated file (tsv) with columns:
                'chrom', 'start', 'end', BAM_FILE_1, BAM_FILE_2 ...
    """
    if type(bams) is str:
        bams = [bams]

    bam_str = ' '.join(bams)
    temp = get_temp_path(prefix='bedtools_multicov_')
    call(f'bedtools multicov -bams {bam_str} -bed {bed} > {temp}')

    # Write header line (columns)
    with open(temp, 'r') as reader:
        with open(output, 'w') as writer:
            bam_columns = '\t'.join(bams)
            writer.write(f'chrom\tstart\tend\t{bam_columns}\n')
            writer.write(reader.read())

    call(f'rm {temp}')
