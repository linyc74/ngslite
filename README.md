[![Coverage Status](https://coveralls.io/repos/github/linyc74/ngslite/badge.svg?branch=master)](https://coveralls.io/github/linyc74/ngslite?branch=master)

# ngslite
**Light-weight functions for next-generation sequencing (NGS) data analysis**

## Install

    pip install ngslite

## Common tools

Common wrapper functions for command-line tools. `samtools` needs to be installed.

    import ngslite as ngs

    ngs.sam_to_bam('path/to/sam')
    ngs.sort_bam('path/to/bam')
    ngs.index_bam('path/to/bam')

    # Or simply combine the three steps above
    ngs.sam_to_indexed_bam('path/to/sam')

Trimming reads. `trim_galore` needs to be installed.

    ngs.trim_galore(fq1='reads.1.fq', fq2='reads.2.fq', gzip=True)
        
Mapping reads. `bowtie2` needs to be installed.

    ngs.bowtie2_mapping(ref='ref.fa', fq1='reads.1.fq', fq2='reads.2.fq', sam='output.sam')

## Handling files

### Fasta

Read fasta files:

    from ngslite import FastaParser

    with FastaParser('file.fa') as parser:
        for header, sequence in parser:
            print(header, sequence)

Read the whole fasta file at once:

    from ngslite import read_fasta
    
    fasta_data = read_fasta('file.fa')
    header, sequence = fasta_data[0]
    print(header, sequence)

Write fasta files:

    from ngslite import FastaWriter

    with FastaWriter('file.fa') as writer:
        writer.write(header, sequence)

### Genbank

Read Genbank files:

    from ngslite import read_genbank

    chromosomes = read_genbank('file.gbk')
    for chromosome in chromosomes:
        print(chromosome.seqname)
        for feature in chromosome.features:
            print(feature)

### GFF

Read GFF files:

    from ngslite import GffParser

    with GffParser('file.gff') as parser:
       for feature in parser:
           print(feature.start,
                 feature.end,
                 feature.strand,
                 feature.attributes)

Write GFF files:

    from ngslite import GffFeature, GffWriter

    feature = GffFeature(
        seqid='chr1',
        source='.',
        type='CDS',
        start=1,
        end=99,
        score='.',
        strand='+',
        phase='.',
        attributes='name=...'
    )

    with GffWriter('file.gff') as writer:
        writer.write(feature)
