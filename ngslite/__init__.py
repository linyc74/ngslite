"""
Light-weight functions for next-generation sequencing (NGS) data analysis

Python API for NGS-related command line tools
Python functions for manipulating NGS-related data and files
"""


__version__ = '0.4'


from .arrfunc import *
from .bedtools import *
from .bowtie2 import *
from .bwa import *
from .cmdtools import *
from .count import *
from .dnatools import *
from .fasta import *
from .fastq import *
from .file_conversion import *
from .filetools import *
from .gtftools import *
from .hmmer import *
from .hmmer_parse import *
from .jellyfish import *
from .kmertools import *
from .metaspades import *
from .samtools import *
from .subset import *
from .trim import *
from .trinity import *
