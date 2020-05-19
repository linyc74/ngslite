"""
Light-weight functions for next-generation sequencing (NGS) and genomics data analysis

Python API for NGS-related command line tools
Python functions for manipulating NGS and genomics data
"""


__version__ = '1.1.1'


from .cmd_toolkit import bedtools_multicov, bowtie2_mapping, bwa_mapping, \
    metaspades, prokka, trim_galore, orf_finder
from .count import count_reads
from .dataclass import GenericFeature, FeatureArray, Chromosome
from .dnatools import rev_comp, translate, base_content
from .fasta import FastaParser, FastaWriter, read_fasta, write_fasta, subset_fasta
from .fasta_gtf import read_fasta_gtf
from .fastq import FastqParser, FastqWriter, interleave
from .file_conversion import sam_to_bam, bam_to_sam, fq_to_fa, vcf_to_bcf, bcf_to_vcf
from .filetools import get_files, change_extension, change_prefix, concat, zip_broadcast
from .genbank_parse import read_genbank, genbank_to_fasta, genbank_to_gtf
from .genbank_write import write_genbank, make_genbank
from .gfftools import GffParser, GffWriter, read_gff, write_gff, subset_gff, print_gff
from .glimmer import glimmer3
from .gtftools import GtfParser, GtfWriter, read_gtf, write_gtf, subset_gtf, print_gtf
from .hmmer import hmmsearch, hmmbuild, parse_hmmsearch_result, validate_hmm_parse_result
from .locus_extractor import locus_extractor, genbank_locus_extractor
from .lowlevel import call, gzip
from .merge_pfam_orf import merge_pfam_into_orf
from .multiple_sequence_alignment import muscle, simple_phylogeny, draw_tree
from .samtools import SamParser, SamWriter, sort_bam, index_bam, sam_to_indexed_bam, \
    subset_bam_regions, remove_unmapped, decode_flag, encode_flag, filter_sam_by_flag, \
    print_flag, print_sam
from .random import random_sample
from .synteny import synteny
from .vcftools import VcfParser, VcfWriter, bcftools_variant_call, sort_bcf, \
    subset_bcf_regions, print_vcf, unpack_vcf_info
