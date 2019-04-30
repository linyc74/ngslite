from .gtftools import read_gtf, write_gtf
from .data_class import FeatureArray, GenericFeature


def _pfam_in_orf(pfam, orf):
    """
    Args:
        pfam: GenericFeature object

        orf: GenericFeature object

    Returns: bool
        Whether the pfam domain is part of the orf
    """
    return pfam.start >= orf.start and \
           pfam.end <= orf.end and \
           pfam.strand == orf.strand and \
           (pfam.start - orf.start) % 3 == 0  # in-frame


def _merge_pfam_arr_to_orf_arr(pfam_arr, orf_arr):
    """
    Args:
        pfam_arr: list of GenericFeature
            Each GenericFeature is a Pfam domain

        orf_arr: list of GenericFeature
            Each GenericFeature is an ORF

    Returns: list of GenericFeature
        Each GenericFeature is an ORF with Pfam information appended
    """
    merge_arr = []

    # Note that GenericFeature objects in orf_arr and pfam_arr will be altered
    for orf in orf_arr:
        while len(pfam_arr) > 0:
            pfam = pfam_arr[0]  # the first pfam feature

            # 1) Pfam is contained by the ORF, so add Pfam information into the ORF
            if _pfam_in_orf(pfam, orf):
                pfam_attr = pfam.attributes[:]
                # Add positional information of the pfam domain
                pfam_attr += [('start', pfam.start), ('end', pfam.end), ('strand', pfam.strand)]
                # Convert all pfam attributes [(key, val)] --> str --> 'note' in the ORF
                orf.add_attribute(key='note', val=str(pfam_attr))

                pfam_name = pfam.get_attribute('name')
                orf_name = orf.get_attribute('name')

                orf.set_attribute('name', f'{orf_name} | {pfam_name}')
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

        merge_arr.append(orf)

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
    pfam_dict = read_gtf(file=pfam, as_dict=True)
    orf_dict = read_gtf(file=orf, as_dict=True)
    merge_dict = {}

    # For each contig, merge Pfam feature array into ORF feature array
    for seqname in orf_dict.keys():

        orf_arr = FeatureArray(
            seqname,
            genome_size=1e6,
            features=orf_dict[seqname],
            circular=False
        )
        orf_arr.sort()

        pfam_arr = FeatureArray(
            seqname,
            genome_size=1e6,
            features=pfam_dict.get(seqname, []),
            circular=False
        )
        pfam_arr.sort()

        merge_dict[seqname] = _merge_pfam_arr_to_orf_arr(pfam_arr, orf_arr)

    write_gtf(data=merge_dict, file=output)
