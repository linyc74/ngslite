from .fasta import FastaParser
from .bedtools import bedtools_multicov
from .lowlevel import __temp
import os
import pandas as pd


def compute_enrichment(contig_fasta, bamC, bamE, totalC, totalE,
                       nameC='control_(RPKM)', nameE='enzyme_(RPKM)', output_csv='enrichment.csv'):
    """
    Use mapped reads to compute enrichment (enzyme / control)

    For each contig, enrichment is calculated as:

                      # mapped enzyme reads (RPKM)
        Enrichment = ------------------------------
                      # mapped control reads (RPKM)

        RPKM = Reads Per Kb per Million reads
               `     `   `      `

    Args:
        contig_fasta: str, path-like
            The fasta containing all contigs

        bamC: str, path-like
            The "Control" bam file which contains reads mapped to the contigs in <contig_fasta>

        bamE: str, path-like
            The "Enzyme" bam file which contains reads mapped to the contigs in <contig_fasta>

        totalC: int
            Total number of reads of "Control"

        totalE: int
            Total number of reads of "Enzyme"

        nameC: str
            Name of the "Control" sample

        nameE: str
            Name of the "Enzyme" sample

        output_csv: str, path-like
            The output csv file containing data of all contigs
    """
    # Generate temp file names
    bed = __temp('temp', '.bed')
    multicov = __temp('multicov', '.tsv')

    # Generate a BED file from all contigs
    with open(bed, 'w') as fh:
        with FastaParser(contig_fasta) as parser:
            for head, seq in parser:
                fh.write(f"{head}\t{0}\t{len(seq)}\n")

    # Count reads mapped to each contig with "bedtools multicov"
    bedtools_multicov(
        bed=bed,
        bams=[bamC, bamE],
        output=multicov
    )

    # Read multicov results as a dataframe
    df = pd.read_csv(multicov, sep='\t')
    df = df.rename(
        columns={'chrom': 'contig_id',
                 bamC: nameC,
                 bamE: nameE}
    )

    # Remove temp files
    os.remove(multicov)
    os.remove(bed)

    # Compute RPKM
    df[nameC] = df[nameC] / (df['end']/1e3) / (totalC/1e6)
    df[nameE] = df[nameE] / (df['end']/1e3) / (totalE/1e6)

    # Compute enrichment
    df['enrichment'] = df[nameE] / df[nameC]
    df.to_csv(output_csv, header=True, index=False)

