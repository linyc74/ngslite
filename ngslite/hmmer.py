from .fasta import FastaParser, FastaWriter, read_fasta
from .dnatools import translate, rev_comp
from .lowlevel import __call
from functools import partial
printf = partial(print, flush=True)


def __gzip(file, keep=True):
    """
    Call the "gzip" command to zip or unzip files.

    Args:
        file: str, path-like

        keep: bool, keep the input file or not
    """
    decomp = ['', '-d '][file.endswith('.gz')]
    keep = ['', '-k '][keep]
    __call(f"gzip {decomp}{keep}{file}")


def __is_dna(fasta):
    """
    Returns 'True' if the first sequence of the <fasta> file is a DNA
    """
    with FastaParser(fasta) as parser:
        is_dna = True
        for char in set(parser.next()[1]):
            if not char in ('A', 'C', 'G', 'T'):
                is_dna = False
                break
    return is_dna


def __translate_dna_database(fasta):
    """
    Translate the input fasta file in six frames.
    Append ';frame=<frame>' in the header line.
    Write a new '<fasta>_translated.fa' file.

    Args:
        fasta: str, path-like

    Returns: str,
        The written fasta name '<fasta>_translated.fa'
    """
    if fasta.endswith('.fa'): output = fasta[:-3] + '_translated.fa'
    elif fasta.endswith('.fasta'): output = fasta[:-6] + '_translated.fa'

    with FastaParser(fasta) as parser:
        with FastaWriter(output, 'w') as writer:
            for head, seq in parser:
                for frame in (1, 2, 3):
                    new_head = f"{head};frame={frame}"
                    aa_seq = translate(seq[frame-1:])
                    writer.write(new_head, aa_seq)

                for frame in (-1, -2, -3):
                    new_head = f"{head};frame={frame}"
                    rc_seq = rev_comp(seq)
                    aa_seq = translate(rc_seq[-frame-1:])
                    writer.write(new_head, aa_seq)

    return output


def __check_fasta_header(fasta):
    """
    Fasta headers should NOT contain ' '. This function checks if it's correct.

    Args:
        fasta: str, path-like

    Returns: bool
    """
    is_correct = True
    with FastaParser(fasta) as parser:
        for head, seq in parser:
            if ' ' in head:
                is_correct = False
                break
    return is_correct


def hmmsearch(hmm, database, output, cpu=2):
    """
    Args:
        hmm: str, path-like
            The .hmm file built by the command "hmmbuild",
                e.g. "hmmbuild -o Pfam-A_summary.txt Pfam-A.hmm Pfam-A.seed"
            Also accepts .gz format

        database: str, path-like
            The fasta database to be searched against
            Also accepts .gz format

        output: str, path-like
            The output text file reported by the command "hmmsearch"

        cpu: int,
            Number of CPUs, 2 is the default in hmmsearch
    """
    db_is_gz = False
    if database.endswith('.gz'):
        __gzip(database, keep=True)
        database = database[:-len('.gz')]
        db_is_gz = True

    hmm_is_gz = False
    if hmm.endswith('.gz'):
        __gzip(hmm, keep=True)
        hmm = hmm[:-len('.gz')]
        hmm_is_gz = True

    db_is_dna = False
    if __is_dna(database):
        database = __translate_dna_database(database)
        db_is_dna = True

    if not __check_fasta_header(database):
        printf(f"Because the database '{database}' contains blank space in headers, hmmsearch is aborted")
        return

    # Run hmmsearch
    __call(f"hmmsearch --cpu {cpu} {hmm} {database} > {output}")

    if db_is_gz or db_is_dna:
        __call(f"rm {database}")
    if hmm_is_gz:
        __call(f"rm {hmm}")


def hmmbuild(seed, hmm):
    """
    Args:
        seed: str, path-like
            The input .seed file

        hmm: str, path-like
            The output .hmm file
    """
    if not hmm.endswith('.hmm'):
        hmm = hmm + '.hmm'
    summary = hmm[:-len('.hmm')] + '_summary.txt'

    __call(f"hmmbuild -o {summary} {hmm} {seed}")


def parse_hmmsearch_result(file, output, database=None):
    """
    Parse the result reported by HMMER, to create a GTF file.

    The result text file output by hmmsearch is highly unstructured.

    For each hmm query:
        Get <query_name>, <query_accession>, <query_description>

        For each translated contig sequence:
            Get <contig_id>, <contig_length_bp>, <frame>

            For each hit (of the hmm domain) with satisfied threshold:
                1) Get <start_aa> and <end_aa> of amino acid position
                2) Use <start_aa>, <end_aa>, <frame> and <contig_length_bp> to compute <start_bp>, <end_bp>

                Write the domain into GTF as a line, consisting of 9 tab-delimited fields:
                1   seqname     <contig_id>
                2   source      .
                3   feature     CDS
                4   start       <start_bp>
                5   end         <end_bp>
                6   score       .
                7   strand      +/-
                8   frame       0
                9   attribute   name "<query_name>";accession "<query_accession>";description "<query_description>";E_value "<E_value>"

    Args:
        file: str, path-like
            The text file reported by "hmmsearch"

        output: str, path-like
            The output GTF file

        database: str, path-like
            The fasta database used for hmmsearch
            This is used to get the contig length.
            If None, then use the contig name in the input <file> to get contig length
    """
    # Read the hmmsearch result
    with open(file, 'r') as fh:
        text = fh.read()

    # Set up a dictionary of contig lengths
    if database:
        with FastaParser(database) as parser:
            contig_length_dict = {head: len(seq) for head, seq in parser}
    else:
        contig_length_dict = {}

    # The output GTF file
    gtf = open(output, 'w')
    gtf_line_count = 0

    # Remove the document header
    text = text[text.find('\n\n') + 2:]

    # Use '//' to split into queries, each query is the result from a particular Pfam-A domain
    all_queries = text.split('\n//\n')[:-1]

    # Remove the sections without hits
    queries = []
    for query in all_queries:
        if not '[No hits detected that satisfy reporting thresholds]' in query:
            queries.append(query)

    # First layer of for loop: Each hmm query is a section
    for query in queries:
        # Get the first three lines
        line1, line2, line3 = query.splitlines()[:3]

        ### Three variables for the query ###
        query_name = line1.split(':')[1].strip()
        query_accession = line2.split(':')[1].strip()
        query_description = line3.split(':')[1].strip()
        ### Three variables for the query ###

        # Get the text between 'Domain annotation for...' and 'Internal pipeline...'
        text = query.split('Domain annotation for each sequence (and alignments):\n')[1].split(
            'Internal pipeline statistics summary:')[0]

        # Use '>> ' to split into translated contigs
        translated_contigs = text.split('>> ')[1:]

        # Second layer of for loop: Each translated contig hit by the query
        for contig in translated_contigs:
            # Get the table appearing before 'Alignments for each domain'
            table = contig.split('\n\n  Alignments for each domain:')[0]

            # line1 is the header line for the translated contig
            line1 = table.splitlines()[0]
            contig_id, frame = line1.split(';frame=')
            frame = int(frame)
            if contig_length_dict:
                contig_length_bp = contig_length_dict[contig_id]
            else:
                contig_length_bp = int(contig_id.split('len=')[1])

            # From the 4th line, each line is a hit
            hits = table.splitlines()[3:]

            # Third layer of for loop: each hit of a query within the translated contig
            for hit in hits:
                number, satisfied, score, \
                bias, c_Evalue, i_Evalue, \
                hmm_from, hmm_to, symbol_1, \
                ali_from, ali_to, symbol_2, \
                env_from, env_to, symbol_3, acc = hit.split()

                # '!' means: The hit satisfies both per-sequence (contig) and per-domain (query) inclusion thresholds
                if satisfied == '!':
                    # Both <start_aa> and <end_aa> are 1-based index
                    start_aa = int(ali_from)  # aligned from
                    end_aa = int(ali_to)  # aligned to
                    i_Evalue = float(i_Evalue)  # independent E-value

                    # All the bp index are 1-based for GTF files
                    # If frame = 1 or -1, then add (abs(frame) - 1) = 0
                    # Do NOT forget to take absolute value of frame! I spent hours on this bug.
                    start_bp = ((start_aa - 1) * 3 + 1) + (abs(frame) - 1)
                    end_bp = (end_aa * 3) + (abs(frame) - 1)

                    # If + strand, then do nothing to the start_bp and end_bp
                    if frame > 0:
                        strand = '+'
                    # If - strand, then invert it from the - strand positions to the + strand positions
                    else:
                        strand = '-'
                        start_bp, end_bp = contig_length_bp - end_bp + 1, contig_length_bp - start_bp + 1

                    # Get attributes for the particular query, i.e. the Pfam domain
                    attribute = 'name "{}";accession "{}";description "{}";E_value "{}"'.format(
                        query_name, query_accession, query_description, i_Evalue)

                    line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        contig_id,  # seqname
                        '.',        # source
                        'CDS',      # feature
                        start_bp,   # start
                        end_bp,     # end
                        '.',        # score
                        strand,     # strand
                        0,          # frame
                        attribute)  # attribute

                    gtf.write(line)
                    gtf_line_count += 1

    gtf.close()
    printf('There are totally {} hits exported into the GTF file "{}".'.format(gtf_line_count, output))


def validate_hmm_parse_result(gtf, dna_database, output):
    """
    Parsing the HMMER result was very complicated.
    In particular, the amino acid positions were converted back to the nucleotide positions.
    There could be errors of this conversion. So I want to verify the nucleotide positions of each predicted domain is correct.

    To verify the positions, for each hit in the GTF file:
        Get DNA sequences from the original fasta database
        Translate it into proteins
        See if there's any stop codon *

    After running the script, I saw there are still a few translated sequences that have stop codon.
    I went into the HMMER result file and found there was indeed stop codon within the predicted domain.

    Args:
        gtf: str, path-like
            The GTF file exported by parse_hmmsearch_result()

        dna_database: str, path-like
            The original DNA database (fasta) used from hmmsearch

        output: str, path-like
            The output text file summarizing the validation result
    """
    contigs = read_fasta(dna_database)

    # Build a dictionary of contigs {header: sequence, ...}
    contigs = {head: seq for head, seq in contigs}

    count_plus = 0
    count_plus_err = 0
    count_minus = 0
    count_minus_err = 0
    with open(gtf, 'r') as gtf_fh:
        for line in gtf_fh:
            seqname, source, feature, start, end, score, strand, frame, attribute = \
                line.rstrip().split('\t')
            start, end = int(start), int(end)

            cds = contigs[seqname][start - 1:end]

            if strand == '+':
                count_plus += 1
            else:
                cds = rev_comp(cds)
                count_minus += 1

            if '*' in translate(cds):
                if strand == '+':
                    count_plus_err += 1
                else:
                    count_minus_err += 1

    with open(output, 'w') as fh:
        text = f"""\
From the GTF file {gtf}:

    There are totally {count_plus} (+) strand CDS in the database {dna_database}, of which:
        {count_plus_err} contains stop codon (*) in the translated sequence.

    There are totally {count_minus} (-) strand CDS in the database {dna_database}, of which:
        {count_minus_err} contains stop codon (*) in the translated sequence."""

        fh.write(text)

