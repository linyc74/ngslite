class FastqParser:
    """
    A simple fastq parser that parses each read (four lines) of a fastq file
    """
    def __init__(self, file):
        """
        Args:
            file: str, path-like object
        """
        self.__fastq = open(file, 'r')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__fastq.seek(0)
        return self

    def __next__(self):
        r = self.next()
        if r:
             return r
        else:  # r is None:
            raise StopIteration

    def next(self):
        """
        Returns: tuple
            The next read of the fastq file.
            If it reaches the end of the file, return (None,)*4.
        """
        line1 = self.__fastq.readline().rstrip()
        line2 = self.__fastq.readline().rstrip()
        line3 = self.__fastq.readline().rstrip()
        line4 = self.__fastq.readline().rstrip()
        if line1:
            return line1, line2, line3, line4
        else:  # line1 == ''
            return None

    def close(self):
        self.__fastq.close()


class FastqWriter:
    """
    A simple fastq writer that writes a single read (four lines) each time
    """
    def __init__(self, file, mode='w'):
        """
        Args
            file: str, path-like object
            mode: str, 'w' for write or 'a' for append
        """
        self.__fastq = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def write(self, read):
        """
        Args:
            read: tuple of four lines (strings) of each read of fastq
        """
        for line in read:
            self.__fastq.write(line+'\n')

    def close(self):
        self.__fastq.close()


def interleave(fq1, fq2, fq_out):
    """
    Interleave two paired-end fastq files.
    This method does not check the header to see
        if the two fastq files are really from paired end reads.
    It just interleave 4 lines by 4 lines.

    Args:
        fq1: str, path to fastq1 (.1.fq)
        fq2: str, path to fastq2 (.2.fq)
        fq_out:
    """
    fq1 = open(fq1, 'r')
    fq2 = open(fq2, 'r')
    fq_out = open(fq_out, 'w')

    while True:
        read1 = [fq1.readline().rstrip() for _ in range(4)] # 4 lines -> list
        read2 = [fq2.readline().rstrip() for _ in range(4)] # 4 lines -> list

        # The header line should start with '@'
        if not read1[0].startswith('@') or not read2[0].startswith('@'):
            break

        for line in read1:
            fq_out.write(line+'\n')
        for line in read2:
            fq_out.write(line+'\n')

    fq1.close()
    fq2.close()
    fq_out.close()

