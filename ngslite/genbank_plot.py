import matplotlib.pyplot as plt
from typing import List, Optional, Dict
from matplotlib.pyplot import Figure
from dna_features_viewer import GraphicFeature, GraphicRecord
from .genbank_parse import read_genbank
from .dataclass import Chromosome, GenericFeature


COLOR_DICT = {
    'CDS': '#993366',
    'tRNA': '#228b22',
}
DEFAULT_COLOR = '#ffffff'


INCH_PER_KB = 2 / 2.54
INCH_PER_RECORD = 10. / 2.54
THICKNESS = 10.
LINEWIDTH = 0.5
LINECOLOR = '#000000'
FONTDICT = None
OPEN_LEFT = False
OPEN_RIGHT = False
BOX_LINEWIDTH = 1
BOX_COLOR = 'auto'


class PlotGenbank:

    file: str
    output: str
    label_attribute: Optional[str]
    skip_types: List[str]
    figure_width: Optional[float]
    figure_height: Optional[float]
    dpi: int

    color_dict: Dict[str, str] = COLOR_DICT
    default_color: str = DEFAULT_COLOR

    inch_per_kb: float = INCH_PER_KB
    inch_per_record: float = INCH_PER_RECORD
    thickness: float = THICKNESS
    linewidth: float = LINEWIDTH
    linecolor: str = LINECOLOR
    fontdict: Dict[str, str] = FONTDICT
    open_left: bool = OPEN_LEFT
    open_right: bool = OPEN_RIGHT
    box_linewidth: float = BOX_LINEWIDTH
    box_color: str = BOX_COLOR

    def filter(
            self,
            generic_features: List[GenericFeature]) \
            -> List[GenericFeature]:
        return [f for f in generic_features if f.type not in self.skip_types]

    def to_graphic_feature(
            self,
            generic_feature: GenericFeature) -> GraphicFeature:

        f = generic_feature
        strand = +1 if f.strand == '+' else -1

        if self.label_attribute is None:
            label = None
        else:
            label = str(f.get_attribute(self.label_attribute))

        color = self.color_dict.get(f.type, self.default_color)

        return GraphicFeature(
            start=f.start,
            end=f.end,
            strand=strand,
            label=label,
            color=color,
            thickness=self.thickness,
            linewidth=self.linewidth,
            linecolor=self.linecolor,
            fontdict=self.fontdict,
            open_left=self.open_left,
            open_right=self.open_right,
            box_linewidth=self.box_linewidth,
            box_color=self.box_color)

    def to_graphic_record(
            self,
            chromosome: Chromosome) -> GraphicRecord:

        generic_features = list(chromosome.features)

        generic_features = self.filter(generic_features=generic_features)

        graphic_features = list(
            map(self.to_graphic_feature, generic_features))

        return GraphicRecord(
            sequence_length=len(chromosome.sequence),
            sequence=chromosome.sequence,
            features=graphic_features)

    def get_max_seq_length(
            self,
            graphic_records: List[GraphicRecord]) -> int:

        seq_lengths = [len(g.sequence) for g in graphic_records]

        return max(seq_lengths)

    def get_empty_figure(
            self,
            graphic_records: List[GraphicRecord]) -> Figure:

        max_seq_length = self.get_max_seq_length(
            graphic_records=graphic_records)

        if self.figure_width is None:
            width = self.inch_per_kb * max_seq_length / 1000
        else:
            width = self.figure_width

        if self.figure_height is None:
            height = self.inch_per_record * len(graphic_records)
        else:
            height = self.figure_height

        figure, axs = plt.subplots(
            nrows=len(graphic_records),
            ncols=1,
            figsize=(width, height)
        )

        return figure

    def plot_graphic_records(
            self,
            graphic_records: List[GraphicRecord]) -> Figure:

        figure = self.get_empty_figure(graphic_records=graphic_records)
        axes = figure.axes

        for i, record in enumerate(graphic_records):

            is_last = i == len(graphic_records) - 1
            with_ruler = True if is_last else False

            record.plot(
                ax=axes[i],
                draw_line=False,
                with_ruler=with_ruler,
                plot_sequence=False,
                annotate_inline=True
            )

        return figure

    def save(self, figure: Figure):
        figure.savefig(self.output, dpi=self.dpi, bbox_inches='tight')

    def main(
            self,
            file: str,
            output: str,
            label_attribute: Optional[str],
            skip_types: List[str],
            figure_width: Optional[float],
            figure_height: Optional[float],
            dpi: int):

        self.file = file
        self.output = output
        self.label_attribute = label_attribute
        self.skip_types = skip_types
        self.figure_width = figure_width
        self.figure_height = figure_height
        self.dpi = dpi

        chromosomes = read_genbank(file=file)

        graphic_records = list(map(self.to_graphic_record, chromosomes))

        figure = self.plot_graphic_records(graphic_records=graphic_records)

        self.save(figure=figure)


def plot_genbank(
        file: str,
        output: str,
        label_attribute: Optional[str] = 'product',
        skip_types: List[str] = ('gene', 'source'),
        figure_width: Optional[float] = None,
        figure_height: Optional[float] = None,
        dpi: int = 600):

    PlotGenbank().main(
        file=file,
        output=output,
        label_attribute=label_attribute,
        skip_types=skip_types,
        figure_width=figure_width,
        figure_height=figure_height,
        dpi=dpi)
