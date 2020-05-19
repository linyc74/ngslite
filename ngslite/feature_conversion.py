from .dataclass import GenericFeature, GtfFeature, GffFeature


def gtf_to_generic_feature(gtf_feature: GtfFeature) -> GenericFeature:

    assert type(gtf_feature) is GtfFeature

    f = gtf_feature

    attr_list = []
    for a in f.attribute.split(';'):
        key = a.split(' "')[0]
        val = a[len(key)+2:-1]
        attr_list.append((key, val))

    return GenericFeature(
        seqname=f.seqname,
        type_=f.feature,
        start=f.start,
        end=f.end,
        strand=f.strand,
        attributes=attr_list,
        frame=f.frame + 1
    )


def gff_to_generic_feature(gff_feature: GffFeature) -> GenericFeature:

    assert type(gff_feature) is GffFeature

    f = gff_feature
    items = [item for item in f.attributes.split(';') if item]

    attr_list = []
    for item in items:
        key, val = item.split('=')
        attr_list.append((key, val))

    return GenericFeature(
        seqname=f.seqid,
        type_=f.type,
        start=f.start,
        end=f.end,
        strand=f.strand,
        attributes=attr_list,
        frame=1 if f.phase == '.' else f.phase + 1
    )


def generic_to_gtf_feature(generic_feature: GenericFeature) -> GtfFeature:
    """
    Covert GenericFeature to GtfFeature (namedtuple)

    Args:
        generic_feature: GenericFeature object
    """
    assert type(generic_feature) is GenericFeature

    f = generic_feature

    # Pack attributes into a single line of str
    attr_str = ''
    for key, val in f.attributes:
        if type(val) is int or type(val) is float:
            attr_str += f"{key} {val};"
        else:  # type(val) is str -> Add quote ""
            val = val.replace(';', '<semicolon>')  # semicolon is not allowed in the GTF attribute field
            attr_str += f"{key} \"{val}\";"

    return GtfFeature(
        seqname=f.seqname,
        source='.',
        feature=f.type,
        start=f.start,
        end=f.end,
        score='.',
        strand=f.strand,
        frame=f.frame - 1,  # GTF frame is 0, 1, 2
        attribute=attr_str[:-1]  # Remove trailing ';'
    )


def generic_to_gff_feature(generic_feature: GenericFeature) -> GffFeature:
    """
    Covert GenericFeature to GffFeature (namedtuple)

    Args:
        generic_feature: GenericFeature object
    """
    assert type(generic_feature) is GenericFeature

    f = generic_feature

    # Pack attributes into a single line of str
    attr_str = ''
    for key, val in f.attributes:
        # Escape characters not allowed
        if type(val) is str:
            val = val.replace(';', '%3B').replace('=', '%3D')
        attr_str += f"{key}={val};"

    return GffFeature(
        seqid=f.seqname,
        source='.',
        type=f.type,
        start=f.start,
        end=f.end,
        score='.',
        strand=f.strand,
        phase=f.frame - 1,  # GFF frame is 0, 1, 2
        attributes=attr_str[:-1]  # Remove trailing ';'
    )
