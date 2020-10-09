import shutil
from ngslite.genbank_plot import plot_genbank
from .setup import setup_dirs, TestCase


class TestGenbankPlot(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_plot_genbank(self):
        plot_genbank(
            file=f'{self.indir}/input.gbk',
            output=f'{self.outdir}/ouptut.pdf',
            label_attribute=None,
            skip_types=['gene', 'source'],
            figure_width=15,
            figure_height=6,
            dpi=600)
