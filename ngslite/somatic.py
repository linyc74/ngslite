import pandas as pd
from typing import Union, List, Dict, Any, Optional


class SortByChromosomePosition:

    df: pd.DataFrame
    chromosome_column: str
    start_column: str
    end_column: str

    chromosome_to_number: Dict[Any, int]

    def main(
            self,
            df: pd.DataFrame,
            chromosome_column: str,
            start_column: str,
            end_column: Optional[str] = None) -> pd.DataFrame:

        self.df = df
        self.chromosome_column = chromosome_column
        self.start_column = start_column
        self.end_column = end_column

        self.set_chromosome_to_number()
        self.sort()

        return self.df

    def set_chromosome_to_number(self):
        chromosomes = self.df[self.chromosome_column].unique()

        self.chromosome_to_number = {}

        # --- canonical chromosomes --- #
        # the number of canonical chromosomes is not the same for different organisms
        # so dynamically build the dictionary
        for chromosome in chromosomes:

            c = str(chromosome)

            if c.startswith('chr'):
                c = c[3:]  # chrX -> X

            if c.isdigit():
                self.chromosome_to_number[chromosome] = int(c)
            elif c == 'X':
                self.chromosome_to_number[chromosome] = 101
            elif c == 'Y':
                self.chromosome_to_number[chromosome] = 102
            elif c == 'MT':
                self.chromosome_to_number[chromosome] = 103

        # --- non-canonical chromosomes --- #
        non_canonical = [
            c for c in chromosomes
            if c not in self.chromosome_to_number.keys()
        ]
        a = 201
        for chromosome in sorted(non_canonical):
            self.chromosome_to_number[chromosome] = a
            a += 1

    def sort(self):
        self.df['order'] = self.df[self.chromosome_column].map(self.chromosome_to_number)

        by = ['order', self.start_column, self.end_column]
        by = [c for c in by if c is not None]
        ascending = [True] * len(by)

        self.df = self.df.sort_values(
            by=by,
            ascending=ascending
        ).drop(
            columns='order'
        )


class FilterVariantsByConsequence:

    PROTEIN_ALTERING = [
        'frameshift_variant',
        'stop_gained',
        'missense_variant',
        'inframe_deletion',
        'inframe_insertion',
        'stop_lost',
        'start_lost',
        'stop_retained_variant',
        'start_retained_variant',
        'incomplete_terminal_codon_variant',
        'coding_sequence_variant',
        'protein_altering_variant',
    ]
    SPLICING = [
        'splice_polypyrimidine_tract_variant',
        'splice_donor_variant',
        'splice_acceptor_variant',
        'splice_donor_region_variant',
        'splice_donor_5th_base_variant',
        'splice_region_variant',
    ]
    INTRON_INTERGENIC = [
        'intron_variant',
        'intergenic_variant',
    ]
    SILENT = [
        'synonymous_variant',
    ]
    UTR_REGULATORY = [
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',

        'upstream_gene_variant',
        'downstream_gene_variant',
        'regulatory_region_ablation',
        'regulatory_region_variant',
        'TF_binding_site_variant',
        'TFBS_ablation',
    ]
    RNA = [
        'mature_miRNA_variant',
        'NMD_transcript_variant',
        'non_coding_transcript_exon_variant',
        'non_coding_transcript_variant',
    ]
    KEEP_CONSEQUENCES = PROTEIN_ALTERING + SPLICING + UTR_REGULATORY

    REQUIRED_COLUMNS = [
        'Consequence',
    ]

    df: pd.DataFrame

    def main(self, df: pd.DataFrame) -> pd.DataFrame:
        self.df = df.copy()

        self.assert_required_columns()
        self.df['has_consequence'] = self.df['Consequence'].apply(self.has_consequence)
        self.df = self.df[self.df['has_consequence']]
        self.df.drop(columns='has_consequence', inplace=True)

        return self.df

    def assert_required_columns(self):
        for col in self.REQUIRED_COLUMNS:
            if col not in self.df.columns:
                raise ValueError(f'Missing required column: "{col}"')

    def has_consequence(self, s: Union[str, float]) -> bool:
        # to deal with multi-consequences like 'intron_variant,non_coding_transcript_variant'
        if pd.isna(s):
            return False
        for c in self.KEEP_CONSEQUENCES:
            if c in s:
                return True
        return False


class FilterGenesByPrevalence:

    REQUIRED_COLUMNS = [
        'SYMBOL',
        'Sample ID',
    ]

    df: pd.DataFrame
    n_samples: int

    genes_to_keep: List[str]

    def main(self, df: pd.DataFrame, n_samples: int) -> pd.DataFrame:
        self.df = df.copy()
        self.n_samples = n_samples

        self.assert_required_columns()
        self.set_genes_to_keep()
        self.df = self.df[self.df['SYMBOL'].isin(self.genes_to_keep)]

        return self.df

    def assert_required_columns(self):
        for col in self.REQUIRED_COLUMNS:
            if col not in self.df.columns:
                raise ValueError(f'Missing required column: "{col}"')

    def set_genes_to_keep(self):
        ret = []  # all genes of interest

        df = self.df.copy().dropna(subset=['SYMBOL']).set_index('SYMBOL')

        for gene in df.index.unique():
            sample_ids = df.loc[gene, 'Sample ID']

            if type(sample_ids) is str:
                sample_ids = [sample_ids]
            else:  # pd.Series of multiple sample_ids
                sample_ids = sample_ids.unique()  # unique sample_ids

            if len(sample_ids) >= self.n_samples:
                ret.append(gene)

        ret = sorted(ret)

        msg = f'{len(ret)} genes to keep:\n{", ".join(ret)}'
        print(msg)

        self.genes_to_keep = ret
