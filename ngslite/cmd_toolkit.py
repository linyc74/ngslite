import os
from typing import Optional, Union, List, Dict
from .gtftools import GtfWriter
from .lowlevel import call, _temp, printf
from .fasta import FastaParser, FastaWriter


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


def bowtie2_mapping(
        ref: str,
        fq1: str,
        sam: str,
        fq2: Optional[str] = None):
    """
    Wrapper function of bowtie2 mapping

    Args:
        ref: path-like
            The reference fasta

        fq1: path-like
            The read-1 fastq

        fq2: path-like
            The read-2 fastq. If none, use <fq1> for unpaired mapping

        sam: path-like
            The output SAM file
    """
    # Build the .bt2 index files
    call(f"bowtie2-build {ref} ref > {ref}_bowtie2_build.log")

    log = sam[:-4]+'.log'
    if fq2 is None:
        # Unpaired mapping
        cmd = f"bowtie2 -x ref -U {fq1} -S {sam} 2> {log}"
    else:
        # Paired-end mapping
        cmd = f"bowtie2 -x ref -1 {fq1} -2 {fq2} -S {sam} 2> {log}"
    call(cmd)

    # Remove the .bt2 index files
    call('rm *.bt2')


def bwa_mapping(
        ref: str,
        fq1: str,
        sam: str,
        fq2: Optional[str] = None,
        threads: int = 4,
        score: int = 30):
    """
    Wrapper function of BWA mapping

    Args:
        ref: path-like
            The reference fasta

        fq1: path-like
            The read-1 fastq

        fq2: path-like
            The read-2 fastq. If none, use <fq1> for unpaired mapping

        sam: path-like
            The output SAM file

        threads:
            Number of CPUs

        score:
            Don’t output alignment with score lower than <score>
    """
    # Build index
    call(f"bwa index {ref}")

    if fq2 is None:
        # Unpaired
        cmd = f"bwa mem -t {threads} -T {score} {ref} {fq1} > {sam}"
    else:
        # Paired-end
        cmd = f"bwa mem -t {threads} -T {score} {ref} {fq1} {fq2} > {sam}"
    call(cmd)

    # Remove the index files
    call('rm ref.*')


def metaspades(
        fq1: str,
        fq2: str,
        output: str,
        min_contig_length: int = 1000,
        threads: int = 16,
        memory: int = 250):
    """
    Wrapper function of the command "spades.py --meta".

    Args:
        fq1: path-like
            The read-1 fastq file

        fq2: path-like
            The read-2 fastq file

        output: path-like
            The output directory of metaspades
            The output fasta file containing assembled contigs

        min_contig_length:
            Minimum (inclusive) contig length in bp

        threads:
            Number of CPUs

        memory:
            Number of Gb of RAM
    """
    call(f"spades.py --meta -1 {fq1} -2 {fq2} -o {output} --threads {threads} --memory {memory} > {output}.log")

    printf(f"Retrieve contigs (>= {min_contig_length} bp) from {output}/contigs.fasta -> {output}.fa")

    with FastaParser(f"{output}/contigs.fasta") as parser:
        with FastaWriter(f"{output}.fa") as writer:
            for head, seq in parser:
                if len(seq) >= min_contig_length:
                    writer.write(head, seq)


def _rename_contig_id(
        genbank: str,
        contig_dict: Dict[str, str]):
    """
    Args:
        genbank: path-like

        contig_dict:
            {'contig_1_52671bp': 'assembly=control_contigs;id=NODE_400_length_52671_cov_411.055515;len=52671'}
               (short header)        (original header)
    """
    temp_gb = _temp('temp', '.gb')  # for example, temp000.gb
    with open(genbank) as reader:
        with open(temp_gb, 'w') as writer:
            for line in reader:
                if line.startswith('LOCUS       '):
                    key = line.split()[1]
                    writer.write(
                        line.replace(key, contig_dict[key])
                    )
                else:
                    writer.write(line)
    os.remove(genbank)
    os.rename(temp_gb, genbank)


def prokka(
        fasta: str,
        outdir: str,
        kingdom: str,
        locus_tag: str = 'LOCUS',
        proteins: str = '',
        metagenome: bool = True,
        evalue: float = 1e-6,
        threads: int = 4,
        log: Optional[str] = None,
        keep: bool = True):
    """
    A wrapper function for Prokka 1.12 (rapid bacterial genome annotation)

    This wrapper function takes a fasta and outputs a genbank file with the name <outdir>.gb

    Args
        fasta: path-like
            The input fasta to be annotated

        outdir: path-like
            The output directory and the genbank file name

        kingdom:
            'Archaea', 'Bacteria', 'Mitochondria' or 'Viruses'

        locus_tag

        proteins: path-like
            FASTA or GBK file to use as 1st priority

        metagenome:
            Improve gene predictions for highly fragmented genomes

        evalue:
            Similarity e-value cut-off

        threads:
            Number of CPUs to be used (0 = all)

        log: path-like
            The log file for stderr

        keep:
            Whether or not to keep all data in the <outdir>
    """
    # Create a temp fasta file to shorten the header because prokka does not take long headers
    # Also create a dictionary to index the shortened header
    temp_fa = _temp('temp', '.fa')  # for example, temp000.fa
    contig_dict = dict()
    with FastaParser(fasta) as parser:
        with FastaWriter(temp_fa) as writer:
            i = 0
            for head, seq in parser:
                new_head = f"contig_{i}_{len(seq)}bp"
                writer.write(new_head, seq)
                contig_dict[new_head] = head
                i += 1

    # Prokka command line
    meta = ['', '--metagenome '][metagenome]
    evalue = format(evalue, 'f')
    if proteins:
        proteins = f"--proteins {proteins} "
    if log is not None:
        log = f" 2> {log}"
    else:
        log = f" 2> {outdir}.log"

    cmd = f"prokka --outdir {outdir} --prefix {outdir} --locustag {locus_tag} \
--kingdom {kingdom} --evalue {evalue} --cpus {threads} {proteins}{meta}{temp_fa}{log}"

    call(cmd)
    os.remove(temp_fa)

    # In the outdir folder, write a tsv file to index shortened contig headers to the original headers
    with open(f"{outdir}/contig_index.tsv", 'w') as fh:
        for key, val in contig_dict.items():
            fh.write(f"{key}\t{val}\n")

    # Copy the genbank file out from the output directory
    call(f"cp {outdir}/{outdir}.gbk {outdir}.gb")

    # Rename the contig header back to the original headers using <contig_dict>
    _rename_contig_id(f"{outdir}.gb", contig_dict)

    # Remove files in the output dir if keep == False
    if not keep:
        call(f"rm -r {outdir}")


def trim_galore(
        fq1: str,
        fq2: str,
        quality: int = 20,
        gzip: bool = True,
        length: int = 20,
        log: str = 'trim_galore.log'):
    """
    Wrapper function of "trim_galore"

    Use the default settings of "trim_galore",
        e.g. Illumina P5 and P7 adapter sequences are used to trim reads

    Args:
        fq1: path-like
            The read-1 fastq file

        fq2: path-like
            The read-2 fastq file

        quality:
            phred33 score

        gzip:
            Compress the output fastq files or not

        length:
            The minimal length (bp) of reads to be retained

        log: path-like
            Append the stderr of trim_galore to the <log> file
    """
    gzip = ['', '--gzip '][gzip]
    call(f'trim_galore --paired --quality {quality} --phred33 --fastqc --illumina \
{gzip}--length {length} --max_n 0 --trim-n --retain_unpaired {fq1} {fq2} 2>> {log}')


def orf_finder(
        fasta: str,
        output: str,
        min_length: int = 75):
    """
    Use the command line tool NCBI ORFfinder to find Open Reading Frames
      in the input fasta file, and output a GTF file

    Args:
        fasta: path-like

        output: path-like
            The output GTF file

        min_length:
            Minimal length of ORF. 75 is the default of ORFfinder
    """
    # ORFfinder can only take fasta headers <= 50 characters
    # Create a temporary fasta file with short headers
    #   and a dictionary matching the short and the original headers
    temp_fa = _temp('temp', '.fa')  # for example, temp000.fa
    contig_dict = {}
    with FastaParser(fasta) as parser:
        with FastaWriter(temp_fa) as writer:
            for i, contig in enumerate(parser):
                head, seq = contig
                writer.write(str(i), seq)
                contig_dict[str(i)] = head

    # -outfmt 3: output format = feature table
    # -ml [int]: minimal length
    temp_table = _temp('temp', '.table')  # for example, temp000.table
    call(f"ORFfinder -in {temp_fa} -out {temp_table} -outfmt 3 -ml {min_length}")

    os.remove(temp_fa)

    with GtfWriter(output) as writer:
        with open(temp_table) as fh:
            while True:
                line1 = fh.readline().rstrip()
                line2 = fh.readline().rstrip()
                line3 = fh.readline().rstrip()
                if line1 == '':
                    break

                # line1: Get contig_number and orf_id
                contig_number = line1.split('_')[1].split(':')[0]
                orf_id = line1.split('>Feature lcl|')[1].split('_')[0]

                # Get the original contig header in the input fasta file
                header = contig_dict[contig_number]

                # line2: Get start and end (zero-based)
                start, end = line2.split('\t')[0:2]
                start, end = int(start), int(end)

                if start < end:
                    strand = '+'
                else:
                    strand = '-'
                    start, end = end, start

                # The start and end of GTF is 1-based
                feature = \
                    (header              ,  # 1 seqname   str
                     '.'                 ,  # 2 source    str
                     'CDS'               ,  # 3 feature   str
                     start + 1           ,  # 4 start     int
                     end + 1             ,  # 5 end       int
                     '.'                 ,  # 6 score     float
                     strand              ,  # 7 strand    str ('+', '-')
                     0                   ,  # 8 frame     int (0, 1, 2)
                     f"name \"{orf_id}\"")  # 9 attribute str
                writer.write(feature)
        os.remove(temp_table)
