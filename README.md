# ngslite
**Light-weight functions for next-generation sequencing (NGS) data analysis**

## Install

    pip install ngslite

## Handling common files

### Fasta

Read fasta files:

    from ngslite import FastaParser, read_fasta

    with FastaParser('file.fa') as parser:
        for header, sequence in parser:
            print(header, sequence)

    fasta_data = read_fasta('file.fa')  # Read the whole fasta file at once
    header, sequence = fasta_data[0]
    print(header, sequence)

Write fasta files:

    from ngslite import FastaWriter

    with FastaWriter('file.fa') as writer:
        writer.write(header, sequence)

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

### Genbank

Read Genbank files:

    from ngslite import read_genbank

    chromosomes = read_genbank('file.gbk')
    for chromosome in chromosomes:
        print(chromosome.seqname)
        for feature in chromosome.features:
            print(feature) 
    
Write Genbank files:

    with ...

### SAM

Read SAM files:

    with ...
