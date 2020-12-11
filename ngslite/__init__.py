__version__ = '1.2.1'

from .cmd_wrapper.prokka import prokka
from .cmd_wrapper.muscle import muscle
from .cmd_wrapper.seqtk import fq_to_fa
from .cmd_wrapper.bwa import bwa_mapping
from .cmd_wrapper.glimmer import glimmer3
from .cmd_wrapper.orf_finder import orf_finder
from .cmd_wrapper.metaspades import metaspades
from .cmd_wrapper.bowtie2 import bowtie2_mapping
from .cmd_wrapper.trim_galore import trim_galore
from .cmd_wrapper.clustalw import simple_phylogeny
from .cmd_wrapper.bedtools import bedtools_multicov
from .cmd_wrapper.bcftools import bcftools_variant_call, vcf_to_bcf, bcf_to_vcf, sort_bcf, subset_bcf_regions
from .cmd_wrapper.hmmer import hmmsearch, hmmbuild, parse_hmmsearch_result, validate_hmm_parse_result
from .cmd_wrapper.samtools import sort_bam, index_bam, sam_to_indexed_bam, subset_bam_regions, remove_unmapped, sam_to_bam, bam_to_sam

from .count import count_reads
from .random import random_sample
from .lowlevel import call, check_output
from .dna import rev_comp, translate, base_content
from .genbank_write import write_genbank, make_genbank
from .fastq import FastqParser, FastqWriter, interleave
from .dataclass import GenericFeature, FeatureArray, Chromosome
from .vcf import VcfParser, VcfWriter, print_vcf, unpack_vcf_info
from .gff import GffParser, GffWriter, read_gff, write_gff, subset_gff, print_gff
from .fasta import FastaParser, FastaWriter, read_fasta, write_fasta, subset_fasta
from .sam import SamParser, SamWriter, decode_flag, encode_flag, filter_sam_by_flag, print_flag, print_sam
from .filetools import get_files, get_dirs, concat, zip_broadcast, get_temp_path, tar, untar, tar_list, gzip
from .genbank_parse import read_genbank, genbank_to_fasta, genbank_to_gff, GenbankTextParser, GenbankParser, GenbankText

from .legacy.fasta_gtf import read_fasta_gtf
from .legacy.gtf import GtfParser, GtfWriter, read_gtf, write_gtf, subset_gtf, print_gtf, genbank_to_gtf
