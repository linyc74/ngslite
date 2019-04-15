from .fasta import read_fasta
from .dnatools import rev_comp
from .jellyfish import jellyfish_count
import numpy as np
import random
import subprocess
from functools import partial
printf = partial(print, flush=True)


def read_kmers(jf_fa, min_count=1):
    """
    Args:
        jf_fa: str, path-like
            The fasta exported by jellyfish

        min_count: int
            Kmers >= min_count will be included

    Returns:
        two arrays:
            kmers: list of str
                Each str is the sequence of a kmer

            counts: np.ndarray, dtype=np.uint32
                Count of kmers in the jf_fa file
    """
    kmers = []
    counts = []
    with open(jf_fa, 'r') as fh:
        while True:
            line1 = fh.readline().strip() # header = count
            line2 = fh.readline().strip() # sequence
            if line1 == '' or line2 == '':
                break
            c = int(line1[1:])
            if c >= min_count:
                kmers.append(line2)
                counts.append(c)
    return kmers, np.array(counts, dtype=np.uint32)


def kmer_2D_points(jf_fa_1, jf_fa_2, min_count=1):
    """
    Args:
        jf_fa_1: str, path-like
            The fasta exported by jellyfish

        jf_fa_2: str, path-like
            The fasta exported by jellyfish

        min_count: int
            Kmers >= min_count will be included from the jellyfish fasta

    Returns:
        three arrays:
            kmers: list of str
                Each str is the sequence of a kmer

            counts_1: np.ndarray, dtype=np.uint32
                Count of kmers in the jf_fa_1 file

            counts_2: np.ndarray, dtype=np.uint32
                Count of kmers in the jf_fa_2 file
    """
    # --- Read the jf_fa_1 file into a dictionary --- #

    # D1 = {'AAAA': count,
    #       'AAAC': count}
    D1 = {}
    with open(jf_fa_1, 'r') as fh:
        while True:
            line1 = fh.readline().strip()
            line2 = fh.readline().strip()
            if line1 == '' or line2 == '':
                break
            c = int(line1[1:])
            if c >= min_count:
                D1.setdefault(line2, c)

    # --- Read the jf_fa_2 file into a dictionary --- #

    # D2 = {'AAAA': count,
    #       'AAAC': count}
    D2 = {}
    with open(jf_fa_2, 'r') as fh:
        while True:
            line1 = fh.readline().strip()
            line2 = fh.readline().strip()
            if line1 == '' or line2 == '':
                break
            c = int(line1[1:])
            if c >= min_count:
                D2.setdefault(line2, c)

    # --- Combine the two dictionaries into one --- #

    # D3 = {'AAAA': (count_1, count_2),
    #       'AAAC': (count_1, count_2)}
    D3 = {}
    for kmer, c1 in D1.items():
        c2 = D2.get(kmer, 0)
        D3.setdefault(kmer, (c1, c2))
        # If the kmer exit in D2, then pop it out
        if c2 != 0:
            D2.pop(kmer)

    # Anything left in D2 is NOT present in D1
    for kmer, c2 in D2.items():
        D3.setdefault(kmer, (0, c2))

    kmers, counts_1_2 = zip(*D3.items())
    counts_1, counts_2 = zip(*counts_1_2)

    return kmers, np.array(counts_1, dtype=np.uint32), np.array(counts_2, dtype=np.uint32)


def save_kmer_2D_points(kmers, x, y, file, x_name='x', y_name='y', description=''):
    """
    This method writes a tab-delimited file to store 2D points of kmers.

    The first line is the description text.

    The second line is the header line:
        'k-mer' \t <x_name> \t <y_name>

    From the third line, each line is a kmer:
        <kmer_seq> \t <x_coordinate> \t <y_coordinate>

    Here's an example:
        k-mer  \t <x_name> \t <y_name>
        AAAAAA \t 10       \t 10
        AAAAAC \t 1        \t 1

    Args:
        kmers: list of str
            Each str is the sequence of a kmer

        x: list of num, or numpy array, dtype = np.int or np.float
            Each num is the x coordinate of the kmer in the <kmers> list

        y: list of num, or numpy array, dtype = np.int or np.float
            Each num is the y coordinate of the kmer in the <kmers> list

        file: str, path-like
            The output tab-delimited file

        x_name: str
            The column name for the x coordinates.
            Usually it's the source filename from which k-mers were counted

        y_name: str
            The column name for the y coordinates.
            Usually it's the source filename from which k-mers were counted

        description: str
            Extra description of this file, which will be written in the first line

    Returns: None
    """
    if len(kmers) != len(x):
        printf("The number of 'kmers' is not equal to that of 'x'")
        return
    if len(kmers) != len(y):
        printf("The number of 'kmers' is not equal to that of 'y'")
        return

    if not file.endswith('.tsv'):
        file = file + '.tsv'

    # Replace any '\n' in the <description> with '   '
    #   since <description> should be limited to the first line
    if '\n' in description:
        description = '   '.join(description.split('\n'))

    with open(file, 'w') as fh:
        fh.write(description + '\n')
        fh.write('k-mer\t{}\t{}\n'.format(x_name, y_name))
        for i in range(len(kmers)):
            fh.write('{}\t{}\t{}\n'.format(kmers[i], x[i], y[i]))


def read_kmer_2D_points(file, dtype=np.uint32, fraction=1.):
    """
    This function reads the tab-delimited file written by the method <save_kmer_2D_points>.

    The first line is the comment/description text.

    The second line is the header line:
        'k-mer' \t <x_name> \t <y_name>

    From the third line, each line is a kmer:
        <kmer_seq> \t <x_coordinate> \t <y_coordinate>

    Here's an example:
        k-mer  \t <x_name> \t <y_name>
        AAAAAA \t 10       \t 10
        AAAAAC \t 1        \t 1

    Args:
        file: str, path-like
            The input tsv file written by the method <save_kmer_2D_points>

        dtype: numpy data types
            Number data types for x and y coordinates; np.uint32 is the default

        fraction: float
            The fraction of k-mers to be randomly retrieved from the file
            Default = 1., i.e. all k-mers

    Returns:
        three arrays:
            kmers: list of str
                Each str is the sequence of a kmer

            x: np.ndarray, dtype=np.uint32
                x coordinate of the kmer in the <kmers> list

            y: np.ndarray, dtype=np.uint32
                y coordinate of the kmer in the <kmers> list

        three strings:
            x_name: str
                The column name for the x coordinates.
                Usually it's the source filename from which k-mers were counted.

            y_name: str
                The column name for the y coordinates.
                Usually it's the source filename from which k-mers were counted.

            description: str
                Extra description of this file
    """
    data = []

    with open(file, 'r') as fh:
        description = fh.readline().rstrip()
        _, x_name, y_name = fh.readline().rstrip().split('\t')
        for line in fh:
            if random.uniform(0, 1) <= fraction:
                data.append(
                    tuple(
                        line.rstrip().split('\t') # 'kmer\tx\ty' => (kmer, x, y)
                    )
                )

    kmers, x, y = zip(*data) # Unzip into three lists

    return kmers, np.array(x, dtype=dtype), np.array(y, dtype=dtype), x_name, y_name, description


def filter_kmer_2D_points(kmers, x, y, ref_fasta=None, min_count=None, func=None):
    """
    Retrieve kmers that fits the following criteria:
        1) Appearing in the reference genome <ref_fasta>. If ref genome is None, then no filtering.
        2) Greater or equal than min_count. If min_count is None, then no filtering.
        3) Passing the criteria of the function <func>. If func is None, then no filtering.

    Args:
        kmers: list of str
            Each str is the sequence of a kmer

        x: np.ndarray, dtype = int or float
            x coordinate of the kmer in the <kmers> list

        y: np.ndarray, dtype = int or float
            y coordinate of the kmer in the <kmers> list

        ref_fasta: str
            If str, then it's the fasta file path of the reference genome
            If None, then don't filter

        min_count: single number or a tuple of (number, number), where the number could be int or float
            If a single number, then filter both x and y with that number
            If (number_x, number_y), then filter x and y with number_x and number_y, respectively
            If None, then don't filter

        func: function object,
            This function should have three input parameters, for <kmers>, <x> and <y>, respectively;
                and should return True/False
            If None, then don't filter

    Returns:
        three arrays:
            [1] list of str
                Each str is the sequence of a kmer

            [2] np.ndarray, dtype determined by the input
                x coordinate of the kmer in the <kmers> list

            [3] np.ndarray, dtype determined by the input
                y coordinate of the kmer in the <kmers> list
    """
    # From a high level, the function walks through the input kmer list,
    #     for each kmer see if it fits the criteria
    #     and eventually return the k-mers fitting the criteria

    kmers_new, x_new, y_new = kmers, x, y

    # --- Filter by the reference genome sequence --- #
    if not ref_fasta is None:

        # Push the new ones into old ones and empty the new ones
        kmers_old, x_old, y_old = kmers_new, x_new, y_new
        kmers_new, x_new, y_new = [], [], []

        ref_seq = read_fasta(ref_fasta)[0][1]

        # Convert to UPPER CASE!
        ref_seq = ref_seq.upper()

        # Build a dictionary of k-mers of the reference sequence: {'AAAA': True, 'AAAC': True}
        # A particular k-mer may appear in the reference sequence more than once,
        #     but it would be unique as a key in the dictionary
        ref_dict = {}

        k = len(kmers[0])
        # For each k-mer in the forward sequence of the reference
        for i in range(len(ref_seq) - k + 1):
            ref_dict.setdefault(ref_seq[i:i+k], True)

        # Reverse complementary of the reference genome sequence
        rc_ref_seq = rev_comp(ref_seq)

        # For each k-mer in the reverse complementary sequence of the reference
        for i in range(len(rc_ref_seq) - k + 1):
            ref_dict.setdefault(rc_ref_seq[i:i+k], True)

        # Now the k-mer dictionary of the reference sequence has been built

        # For each input k-mer, see if it appears in the reference genome
        # If yes
        for i, kmr in enumerate(kmers_old):
            if ref_dict.get(kmr, None):
                kmers_new.append(kmr)
                x_new.append(x_old[i])
                y_new.append(y_old[i])

    # --- Filter by minimum count --- #
    if not min_count is None:

        # Push the new ones into old ones and empty the new ones
        kmers_old, x_old, y_old = kmers_new, x_new, y_new
        kmers_new, x_new, y_new = [], [], []

        if isinstance(min_count, int):
            x_min, y_min = min_count, min_count
        if isinstance(min_count, tuple) or isinstance(min_count, list):
            x_min, y_min = min_count

        for i, kmr in enumerate(kmers):
            if x_old[i] >= x_min and y_old[i] >= y_min:
                kmers_new.append(kmr)
                x_new.append(x_old[i])
                y_new.append(y_old[i])

    # --- Filter by the function object --- #
    if not func is None:

        # Push the new ones into old ones and empty the new ones
        kmers_old, x_old, y_old = kmers_new, x_new, y_new
        kmers_new, x_new, y_new = [], [], []

        for i, kmr in enumerate(kmers):
            if func(kmr, x_old[i], y_old[i]):
                kmers_new.append(kmr)
                x_new.append(x_old[i])
                y_new.append(y_old[i])

    return kmers_new, np.array(x_new), np.array(y_new)


def sample_kmer_2D_points(kmers, x, y, fraction):
    """
    Randomly sample kmers

    Args:
        kmers: list of str
            Each str is the sequence of a kmer

        x: np.ndarray, dtype = int or float
            x coordinate of the kmer in the <kmers> list

        y: np.ndarray, dtype = int or float
            y coordinate of the kmer in the <kmers> list

        fraction: float (0..1)
            The fraction of kmers to be sampled

    Returns:
        three arrays:
            [1] list of str
                Each str is the sequence of a kmer

            [2] np.ndarray, dtype determined by the input
                x coordinate of the kmer in the <kmers> list

            [3] np.ndarray, dtype determined by the input
                y coordinate of the kmer in the <kmers> list
    """
    ret_kmer = []
    ret_x = []
    ret_y = []
    for i, kmr in enumerate(kmers):
        if random.uniform(0, 1) <= fraction:
            ret_kmer.append(kmr)
            ret_x.append(x[i])
            ret_y.append(y[i])
    return ret_kmer, np.array(ret_x), np.array(ret_y)


def build_kmer_dict(kmers, x, y):
    """
    Args:
        kmers: list of str
            Each str is the sequence of a kmer

        x: list of num, or numpy array, dtype = np.int or np.float
            Each num is the x coordinate of the kmer in the <kmers> list

        y: list of num, or numpy array, dtype = np.int or np.float
            Each num is the y coordinate of the kmer in the <kmers> list

    Returns: dict
        {'AAAA': (x_1, y_1),
         'AAAC': (x_2, y_2)}
        where the coordinates x_1, y_1, x_2, y_2 could be int or float
        depending on the input <x> and <y>
    """
    return {kmr: xy for kmr, xy in zip(kmers, zip(x, y))}


def get_kmer_2D_points(seq, kmer_dict):
    """
    Convert the input <seq> into k-mers and
        find the 2D positions of those k-mers in the <kmer_dict>

    Args:
        seq: str
            DNA sequence

        kmer_dict: dict
            {'AAAA': (x_1, y_1),
             'AAAC': (x_2, y_2)}
            where the coordinates x_1, y_1, x_2, y_2 could be int or float

    Returns:
        three arrays:
            kmers: list of str
                Each str is the sequence of a kmer

            x: np.ndarray, dtype determined by the input kmer_dict
                x coordinate of the kmer in the <kmers> list

            y: np.ndarray, dtype determined by the input kmer_dict
                y coordinate of the kmer in the <kmers> list
    """
    # Overall, the function does:
    #   1) Walk through each k-mer in the <seq>
    #   2) If the k-mer appears in the <kmer_dict>, put it in the new dictionary
    #      A particular k-mer could appear more than once in the <seq>,
    #          but it would appear only once in the new dictionary
    #   3) Convert the dictionary of k-mers into three lists: kmers, x, y

    # The length of each k-mer
    k = len(tuple(kmer_dict.keys())[0])

    # Convert to UPPER CASE!
    seq = seq.upper()

    # Build a k-mer dictionary for the input <seq>:
    #     {'AAAA': (x, y), 'AAAC': (x, y), ...}
    # A particular k-mer may appear in the input <seq> more than once,
    #     but it would be unique as a key in the dictionary
    seq_kmer_dict = {}

    # For each k-mer in the input sequence,
    #     see if the k-mer is available in the kmer_dict.
    # If available, then get the counts of the k-mer from kmer_dict
    for i in range(len(seq) - k + 1):

        kmr = seq[i:(i+k)]
        rc_kmr = rev_comp(kmr)

        # Either kmr or rc_kmr should be found in the kmer_dict
        #     because of canonical k-mer counting
        # If kmr is found in kmer_dict,
        #     then there is no need to try to find rc_kmr in the kmer_dict
        if kmer_dict.get(kmr, None):
            seq_kmer_dict.setdefault(kmr, kmer_dict[kmr])

        elif kmer_dict.get(rc_kmr, None):
            seq_kmer_dict.setdefault(rc_kmr, kmer_dict[rc_kmr])

    if seq_kmer_dict:
        kmers, xy = zip(*seq_kmer_dict.items())
        x, y = zip(*xy)
    else:
        kmers, x, y = [], [], []

    return kmers, np.array(x), np.array(y)


def kmers_base_content(kmers, base):
    """
    Args:
        kmers: list of str
            An array of k-mers

        base: str
            e.g., 'A' gives the A content; 'GC' gives the GC content; 'ACGT' should give 1.0

    Returns:
        np.ndarray, dtype=np.float32:
            An array of base contents of the input k-mers
    """
    ret = np.zeros((len(kmers), ), dtype=np.float32)
    for b in set(base):
        ret = ret + np.array([kmr.count(b) / len(kmr) for kmr in kmers], dtype=np.float32)
    return ret


def fastq_to_saved_kmer_2D_points(fq1, fq2, k, output, min_count=1, hash_size='100M', threads=4, description=''):
    """
    The full implementation of the pipeline:
        (1) Input two fastq files
        (2) Count k-mers in each fastq file
        (3) Merge (outer join) the k-mer counts of the two files
        (4) Write a tab-separated file (.tsv) of k-mer 2D points

    Args:
        fq1: str, path-like
            The first fastq file

        fq2: str, path-like
            The second fastq file

        k: uint
            The size of the k-mer in bp

        output: str, path-like
            The output tsv file

        min_count: int
            K-mers with count < min_count will be set to zero

        hash_size: str
            The option for the jellyfish_count() function; '100M' = 100 million, '1G' = 1 billion

        threads: uint
            The number of CPU to be used for the jellyfish_count() function

        description: str
            The description to be written into the output tsv file
    """
    # fastq -> fasta reported by jellyfish
    jellyfish_count(file=fq1, k=k, output='1.fa', min_count=min_count, hash_size=hash_size, threads=threads)
    jellyfish_count(file=fq2, k=k, output='2.fa', min_count=min_count, hash_size=hash_size, threads=threads)

    # Merge two fasta reported by jellyfish -> k-mer 2D points
    kmers, counts_1, counts_2 = kmer_2D_points(jf_fa_1='1.fa', jf_fa_2='2.fa', min_count=min_count)

    # Save k-mer 2D points as the <output> tsv file
    save_kmer_2D_points(kmers=kmers, x=counts_1, y=counts_2, file=output, x_name=fq1, y_name=fq2, description=description)

    # Remove the two fasta files exported by jellyfish
    subprocess.check_call('rm 1.fa 2.fa', shell=True)

