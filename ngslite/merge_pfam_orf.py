from .gtftools import *


def __pfam_in_orf(pfam, orf):
    return pfam.start >= orf.start and \
           pfam.end <= orf.end and \
           pfam.strand == orf.strand and \
           (pfam.start - orf.start) % 3 == 0


def __merge_pfam_arr_to_orf_arr(pfam_arr, orf_arr):
    """
    Args:
        pfam_arr: list of GtfFeature
            Each GtfFeature is a Pfam domain

        orf_arr: list of GtfFeature
            Each GtfFeature is an ORF

    Returns: list of GtfFeature
        Each GtfFeature is an ORF with Pfam information appended
    """
    merge_arr = []

    for orf in orf_arr:
        # D1 is the dictionary of orf attribute
        D1 = attribute_str_to_dict(orf.attribute)

        while len(pfam_arr) > 0:
            # Get the first Pfam
            pfam = pfam_arr[0]

            # D2 is the dictionary of Pfam attribute
            D2 = attribute_str_to_dict(pfam.attribute)
            D2['start'] = pfam.start
            D2['end'] = pfam.end

            # 1) Pfam is completely within the ORF, so add Pfam information into the ORF
            if __pfam_in_orf(pfam, orf):
                D1['name'] += ' | ' + D2['name']
                D1.setdefault('note', '')
                D1['note'] += str(D2) + ' | '
                pfam_arr.pop(0)

            # 2) The first Pfam is ahead of ORF, go for the next ORF
            elif pfam.start > orf.end:
                break

            # 3) The first Pfam partially overlaps with the ORF,
            #      or is completely before the ORF,
            #      or is on the opposite strand
            else:
                # This is an orphan Pfam not contained by any ORF, add it into the output array
                merge_arr.append(pfam_arr.pop(0))

        # Remove trailing ' | '
        if D1.get('note', '').endswith(' | '):
            D1['note'] = D1['note'][:-3]

        # Instantiate a new GtfFeature for the re-annotated ORF
        merge_arr.append(
            GtfFeature(
                seqname=orf.seqname,
                source='.',
                feature='CDS',
                start=orf.start,
                end=orf.end,
                score='.',
                strand=orf.strand,
                frame=0,
                attribute=attribute_dict_to_str(D1),  # D1 is the dictionary of orf attribute
            )
        )

    return merge_arr


def merge_pfam_into_orf(pfam, orf, output):
    """
    Args:
        pfam: str, path-like
            The input GTF containing Pfam annotation

        orf: str, path-like
            The input GTF containing ORFs

        output: str, path-like
            The output GTF in which Pfam is merged into ORFs
    """
    pfam_dict = gtf_to_dict(file=pfam)
    orf_dict = gtf_to_dict(file=orf)
    merge_dict = {}

    # For each contig, merge Pfam feature array into ORF feature array
    for seqname in orf_dict.keys():
        orf_arr = orf_dict[seqname]
        orf_arr = sorted(orf_arr, key=lambda x: x.start)
        pfam_arr = pfam_dict.get(seqname, [])
        pfam_arr = sorted(pfam_arr, key=lambda x: x.start)

        merge_dict[seqname] = __merge_pfam_arr_to_orf_arr(pfam_arr, orf_arr)

    dict_to_gtf(dict_=merge_dict, output=output)
