class GtfParser:
    def __init__(self, file):
        """
        Args:
            file: str, path-like object
        """
        self.__gtf = open(file, 'r')

    def next(self):
        """
        Each line of the GTF file has 9 fields````

            Col	Field       Type
        0   1   seqname     string
        1   2   source      string
        2   3   feature     string
        3   4   start       int
        4   5   end         int
        5   6   score       float
        6   7   strand      string ('+', '-')
        7   8   frame       int (0, 1, 2)
        8   9   attribute   string

        Returns: tuple of str or int
            9 fields of a feature (i.e. a line)
        """
        line = self.__gtf.readline().rstrip()
        if line == '':
            return (None, ) * 9
        else:
            fields = line.split('\t')
            for i in (3, 4, 7):
                if fields[i] != '.':
                    fields[i] = int(fields[i])
            if fields[5] != '.':
                fields[5] = float(fields[5])
            return tuple(fields)

    def close(self):
        self.__gtf.close()


class GtfWriter:
    def __init__(self, file, mode='w'):
        """
        Args:
            file: str, path-like object
            mode: str, 'w' or 'a'
        """
        self.__gtf = open(file, mode)

    def write(self, feature):
        """
        Args:
            feature: tuple of str or int
                Containing 9 fields of a line of GTF file
        """
        self.__gtf.write('\t'.join(map(str, feature)) + '\n')

    def close(self):
        self.__gtf.close()


def subset_gtf(file, seqname, output):
    """
    Args:
        file: str, path-like
            The input gtf file

        seqname: str, or list of str
            Each str is a seqname (chromosome name) to be included

        output: str, path-like
            The output gtf file
    """
    if isinstance(seqname, str):
        seqname = [seqname]

    parser = GtfParser(file)
    writer = GtfWriter(output)

    while True:
        item = parser.next()
        if item[0] is None:
            break

        if item[0] in seqname:
            writer.write(item)

    parser.close()
    writer.close()


def gtf_replace_blank_with(file, s, output):
    """
    Replace blank spaces with <s> in the input GTF file.
    Check the fields that are strings: <seqname>, <source>, <feature>, <attribute>

    Args:
        file: str, path-like
            The input GTF file

        s: str
            The string used to replace blank space

        output:
            The output GTF file
    """
    parser = GtfParser(file)
    writer = GtfWriter(output)

    while True:
        seqname, source, feature, start, end, score, strand, frame, attribute = parser.next()
        if seqname is None:
            break

        seqname = s.join(seqname.split(' '))
        source = s.join(source.split(' '))
        feature = s.join(feature.split(' '))
        attribute = s.join(attribute.split(' '))

        writer.write((seqname, source, feature, start, end, score, strand, frame, attribute))

    parser.close()
    writer.close()

