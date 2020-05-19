from dna_features_viewer import GraphicFeature, GraphicRecord
from .dataclass import Chromosome, GenericFeature
from .genbank_parse import read_genbank


COLOR_DICT = {
    'CDS': '#993366',
    'tRNA': '#228b22',
}


def generic_to_graphic_feature(generic_feature, label_attribute):
    f = generic_feature
    strand = +1 if f.strand == '+' else -1
    label = str(f.get_attribute(label_attribute))
    color = COLOR_DICT.get(f.type, '#ffffff')

    return GraphicFeature(
        start=f.start, end=f.end, strand=strand, label=label, color=color,
        thickness=10, linewidth=0.5, linecolor='#000000', fontdict=None,
        html=None, open_left=True, open_right=False, box_linewidth=1,
        box_color='auto'
    )


def filter_generic_features(generic_features, types, in_types):
    ret = []
    for f in generic_features:
        if (f.type in types) == in_types:
            ret.append(f)
    return ret


def chromosome_to_graphic_record(chromosome, label_attribute, skip_types):

    gen_features = filter_generic_features(
        generic_features=chromosome.features, types=skip_types, in_types=False
    )

    def convert(generic_feature):
        return generic_to_graphic_feature(generic_feature, label_attribute)
    graph_features = list(map(convert, gen_features))

    return GraphicRecord(
        sequence_length=len(chromosome.sequence),
        sequence=chromosome.sequence,
        features=graph_features
    )


def plot_graphic_record(graphic_record, output, fmt, figure_width):
    if figure_width is None:
        figure_width = len(graphic_record.sequence) / 1000
    ax, _ = graphic_record.plot(figure_width=figure_width)
    ax.figure.savefig(output, fmt=fmt, bbox_inches='tight')


def plot_genbank(file, output, fmt='.png', label_attribute='product',
                 hide_types=('gene', 'source'), figure_width=None):

    chromosome = read_genbank(file=file)[0]

    graph_record = chromosome_to_graphic_record(
        chromosome=chromosome, label_attribute=label_attribute,
        skip_types=hide_types
    )

    plot_graphic_record(
        graphic_record=graph_record, output=output, fmt=fmt,
        figure_width=figure_width
    )
