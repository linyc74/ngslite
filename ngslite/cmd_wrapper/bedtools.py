from typing import Union, List
from ..lowlevel import call


def bedtools_multicov(
        bed: str,
        bams: Union[str, List[str]],
        output: str):
    """
    Wrapper function of the command "bedtools multicov -bams <bams> -bed <bed> > <output>"
    Adds a header line to the <output> file

    Args:
        bed: path-like
            The bed file, or any other interval file accepted by the bedtools

        bams: path-like or list of paths
            The bam file, or any other mapped read files accepted by the bedtools

        output: path-like
            The output tab-separated file (tsv) with headers according to the input file type
            Currently supports the header of the following file types:
                bed file: chrom \t start \t end
                gtf file: seqname \t source \t feature \t start \t end \t score \t strand \t frame \t attribute
    """
    if isinstance(bams, list):
        bams = ' '.join(bams)

    call(f"bedtools multicov -bams {bams} -bed {bed} > {output}")

    # Read the data of the bedtools output
    with open(output, 'r') as fh:
        data = fh.read()
    call(f"rm {output}")

    # Write the header line to the output file
    with open(output, 'w') as fh:
        # Replace ' ' with '\t'
        bams = '\t'.join(bams.split(' '))
        if bed.endswith('.bed'):
            fh.write(f"chrom\tstart\tend\t{bams}\n")
        if bed.endswith('.gtf'):
            fh.write(f"seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute\t{bams}\n")
        fh.write(data)