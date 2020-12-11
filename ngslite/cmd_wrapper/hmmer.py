import os
from typing import Optional, Dict, List, Tuple, Any
from ..filetools import gzip
from ..lowlevel import call, printf
from ..dna import translate, rev_comp
from ..fasta import FastaParser, FastaWriter, read_fasta


def hmmbuild(seed: str, hmm: str):
    """
    Args:
        seed: path-like
            The input .seed file

        hmm: path-like
            The output .hmm file
    """
    if not hmm.endswith('.hmm'):
        hmm = hmm + '.hmm'
    summary = hmm[:-len('.hmm')] + '_summary.txt'

    call(f"hmmbuild -o {summary} {hmm} {seed}")


class Hmmsearch:

    hmm: str
    database: str
    output: str
    cpu: int

    unzipped_hmm: Optional[str] = None
    unzipped_db: Optional[str] = None
    translated_db: Optional[str] = None

    def __is_dna(self, fasta: str) -> bool:
        dna_chars = ('A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n')

        with FastaParser(fasta) as parser:
            is_dna = True
            for char in set(parser.next()[1]):
                if char not in dna_chars:
                    is_dna = False
                    break
        return is_dna

    def __six_frame_translate(self, fna: str) -> str:

        output = fna + '_translated.fa'

        with FastaParser(fna) as parser:
            with FastaWriter(output, 'w') as writer:
                for head, seq in parser:
                    for frame in (1, 2, 3):
                        new_head = f"{head};frame={frame}"
                        aa_seq = translate(seq[frame - 1:])
                        writer.write(new_head, aa_seq)

                    for frame in (-1, -2, -3):
                        new_head = f"{head};frame={frame}"
                        rc_seq = rev_comp(seq)
                        aa_seq = translate(rc_seq[-frame - 1:])
                        writer.write(new_head, aa_seq)

        return output

    def set_hmm(self, hmm: str):

        if hmm.endswith('.gz'):
            gzip(hmm, keep=True)
            hmm = hmm[:-len('.gz')]
            self.unzipped_hmm = hmm  # capture the intermediate file

        self.hmm = hmm

    def set_database(self, database: str):

        if database.endswith('.gz'):
            gzip(database, keep=True)
            database = database[:-len('.gz')]
            self.unzipped_db = database  # capture the intermediate file

        if self.__is_dna(fasta=database):
            database = self.__six_frame_translate(fna=database)
            self.translated_db = database  # capture the intermediate file

        self.database = database

    def not_correct_header(self) -> bool:

        with FastaParser(self.database) as parser:
            for head, seq in parser:
                if ' ' in head:
                    return True

        return False

    def print_error_message(self):
        printf(f"Because the database '{self.database}' contains blank space in headers, hmmsearch is aborted")

    def execute(self):
        cmd = f'hmmsearch --cpu {self.cpu} {self.hmm} {self.database} > {self.output}'
        call(cmd)

    def clean_up(self):
        if self.unzipped_hmm is not None:
            os.remove(self.unzipped_hmm)

        if self.unzipped_db is not None:
            os.remove(self.unzipped_db)

        if self.translated_db is not None:
            os.remove(self.translated_db)

    def main(
            self,
            hmm: str,
            database: str,
            output: str,
            cpu: int):

        self.cpu = cpu
        self.output = output

        self.set_hmm(hmm=hmm)
        self.set_database(database=database)

        if self.not_correct_header():
            self.print_error_message()
            return

        self.execute()

        self.clean_up()


def hmmsearch(hmm: str, database: str, output: str, cpu: int = 2):
    """
    Args:
        hmm: path-like
            .hmm file, .gz format accepted

        database: path-like
            The fasta database to be searched against, .gz format accepted

        output: path-like
            The output text file reported by the command "hmmsearch"

        cpu:
            Number of CPUs
    """

    Hmmsearch().main(
        hmm=hmm, database=database, output=output, cpu=cpu)


class ParseSixFrameTranslationHmmsearchTxt:

    file: str
    output: str
    database: str

    contig_len_dict: Dict[str, int]

    gtf: Any  # file io object
    gtf_line_count: int

    def __remove_header(self, full_text: str) -> str:
        t = full_text
        return t[t.find('\n\n') + 2:]

    def get_full_text(self) -> str:
        with open(self.file) as fh:
            ret = self.__remove_header(full_text=fh.read())
        return ret

    def set_contig_len_dict(self):

        if self.database is None:
            self.contig_len_dict = {}

        else:
            with FastaParser(self.database) as parser:
                self.contig_len_dict = {
                    head: len(seq) for head, seq in parser}

    def get_query_blocks(self, full_text: str) -> List[str]:

        blocks = full_text.split('\n//\n')[:-1]

        string = '[No hits detected that satisfy reporting thresholds]'
        blocks = [b for b in blocks if string not in b]

        return blocks

    def get_pfam_info(self, query_block: str) -> Tuple[str, str, str]:

        line1, line2, line3 = query_block.splitlines()[:3]

        pfam_name = line1.split(':')[1].strip()
        pfam_accession = line2.split(':')[1].strip()
        pfam_description = line3.split(':')[1].strip()

        return pfam_name, pfam_accession, pfam_description

    def get_contig_blocks(self, query_block: str) -> List[str]:

        x = 'Domain annotation for each sequence (and alignments):\n'
        y = 'Internal pipeline statistics summary:'

        text_between_x_and_y = query_block.split(x)[1].split(y)[0]

        contig_blocks = text_between_x_and_y.split('>> ')[1:]

        message = '[No individual domains that satisfy reporting thresholds (although complete target did)]'
        contig_blocks = [c for c in contig_blocks if message not in c]

        return contig_blocks

    def get_table_block(self, contig_block: str) -> str:
        return contig_block.split('\n\n  Alignments for each domain:')[0]

    def __get_contig_length(self, contig_id: str) -> int:
        if self.contig_len_dict:
            return self.contig_len_dict[contig_id]
        else:
            return int(contig_id.split('len=')[1])

    def get_contig_info(self, table_block: str) -> Tuple[str, int, int, str]:

        title = table_block.splitlines()[0]

        contig_id, signed_frame = title.split(';frame=')
        signed_frame = int(signed_frame)

        contig_length = self.__get_contig_length(contig_id=contig_id)
        frame = abs(signed_frame)
        strand = '+' if signed_frame > 0 else '-'

        return contig_id, contig_length, frame, strand

    def get_hit_lines(self, table_block: str) -> List[str]:

        hit_lines = []
        for line in table_block.splitlines()[3:]:  # from the 4th line, each line is a hit
            symbol = line.split()[1]

            # Satisfies both per-sequence (contig) and per-domain (query) inclusion thresholds
            if symbol == '!':
                hit_lines.append(line)

        return hit_lines

    def get_hit_info(self, hit_line: str) -> Tuple[int, int, float]:
        items = hit_line.split()

        start_aa = int(items[9])
        end_aa = int(items[10])
        independent_evalue = float(items[5])

        return start_aa, end_aa, independent_evalue

    def __get_reverse_strand_positions(
            self,
            contig_length: int,
            start: int,
            end: int) -> Tuple[int, int]:

        start, end = contig_length - end + 1, contig_length - start + 1

        return start, end

    def aa_to_bp_positions(self, contig_length: int, strand: str, frame: int, start_aa: int, end_aa: int) -> Tuple[int, int]:

        # start_aa: 1-based, inclusive
        # end_bp: 1-based, inclusive
        # bp index for GTF files are 1-based, inclusive
        start_bp = ((start_aa - 1) * 3 + 1) + (frame - 1)
        end_bp = (end_aa * 3) + (frame - 1)

        if strand == '-':
            start_bp, end_bp = self.__get_reverse_strand_positions(
                contig_length=contig_length,
                start=start_bp,
                end=end_bp)

        return start_bp, end_bp

    def get_gtf_attribute(
            self,
            name: str,
            accession: str,
            description: str,
            evalue: float) -> str:

        attribute = f'name "{name}";accession "{accession}";description "{description}";E_value "{evalue}"'

        return attribute

    def write_gtf_line(
            self,
            seqname: str,
            start: int,
            end: int,
            strand: str,
            attribute: str):

        source = '.'
        feature = 'CDS'
        score = '.'
        frame = 0

        items = [
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            frame,
            attribute,
        ]

        line = '\t'.join(map(str, items)) + '\n'

        self.gtf.write(line)
        self.gtf_line_count += 1

    def log(self):
        printf(f'There are totally {self.gtf_line_count} hits exported into the GTF file "{self.output}".')

    def main(
            self,
            file: str,
            output: str,
            database: Optional[str]):

        self.file = file
        self.output = output
        self.database = database

        full_text = self.get_full_text()

        self.set_contig_len_dict()

        self.gtf = open(output, 'w')
        self.gtf_line_count = 0

        query_blocks = self.get_query_blocks(full_text=full_text)  # each query = Pfam domain

        for query_block in query_blocks:

            pfam_name, pfam_accession, pfam_description = self.get_pfam_info(
                query_block=query_block)

            contig_blocks = self.get_contig_blocks(
                query_block=query_block)

            for contig_block in contig_blocks:

                table_block = self.get_table_block(
                    contig_block=contig_block)

                contig_id, contig_length, frame, strand = self.get_contig_info(
                    table_block=table_block)

                hit_lines = self.get_hit_lines(table_block=table_block)

                for hit_line in hit_lines:

                    start_aa, end_aa, independent_evalue = self.get_hit_info(
                        hit_line=hit_line)

                    start_bp, end_bp = self.aa_to_bp_positions(
                        contig_length=contig_length,
                        strand=strand,
                        frame=frame,
                        start_aa=start_aa,
                        end_aa=end_aa)

                    attribute = self.get_gtf_attribute(
                        name=pfam_name,
                        accession=pfam_accession,
                        description=pfam_description,
                        evalue=independent_evalue)

                    self.write_gtf_line(
                        seqname=contig_id,
                        start=start_bp,
                        end=end_bp,
                        strand=strand,
                        attribute=attribute)

        self.gtf.close()

        self.log()


def parse_six_frame_translation_hmmsearch_result(
        file: str, output: str, database: Optional[str] = None):

    ParseSixFrameTranslationHmmsearchTxt().main(
        file=file, output=output, database=database)


def parse_hmmsearch_result(
        file: str, output: str, database: Optional[str] = None):
    """
    Parse the result reported by HMMER, to create a GTF file

    The result text file output by hmmsearch is highly unstructured

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
        file: path-like
            The text file reported by "hmmsearch"

        output: path-like
            The output GTF file

        database: path-like
            The fasta database used for hmmsearch
            This is used to get the contig length
            If None, then use the contig name in the input <file> to get contig length
    """

    parse_six_frame_translation_hmmsearch_result(
        file=file, output=output, database=database)


def validate_hmm_parse_result(
        gtf: str, dna_database: str, output: str):
    """
    Parsing the HMMER result was very complicated
    In particular, the amino acid positions were converted back to the nucleotide positions
    There could be errors of this conversion
    So I want to verify the nucleotide positions of each predicted domain is correct

    To verify the positions, for each hit in the GTF file:
        Get DNA sequences from the original fasta database
        Translate it into proteins
        See if there's any stop codon *

    After running the script, I saw there are still a few translated sequences that have stop codon
    I went into the HMMER result file and found there was indeed stop codon within the predicted domain

    Args:
        gtf: path-like
            The GTF file exported by parse_hmmsearch_result()

        dna_database: path-like
            The original DNA database (fasta) used from hmmsearch

        output: path-like
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
