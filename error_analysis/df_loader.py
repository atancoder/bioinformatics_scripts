import csv
from abc import ABC, abstractmethod

import pandas as pd

from schema import DFSchema


class DFLoader(ABC):
    def __init__(self, filename_tsv: str, tss_filename_bed: str) -> None:
        self.df = pd.read_csv(filename_tsv, delimiter="\t", compression="gzip")
        self.tss_df = self.read_bed_file(tss_filename_bed)

    @classmethod
    def read_bed_file(cls, filename: str) -> pd.DataFrame:
        return pd.read_csv(
            filename,
            sep="\t",
            header=None,
            names=[
                DFSchema.CHROM,
                DFSchema.START,
                DFSchema.END,
                DFSchema.TARGET_GENE,
            ],
            usecols=range(4),
        )

    @abstractmethod
    def load(self) -> pd.DataFrame:
        raise NotImplementedError()

    def _standardize_gene_tss(self, orig_df: pd.DataFrame) -> pd.DataFrame:
        # filter out rows with genes not in tss_df
        orig_df = orig_df[
            orig_df[DFSchema.TARGET_GENE].isin(self.tss_df[DFSchema.TARGET_GENE])
        ].reset_index(drop=True)

        # Use self.tss_df for the target_gene_tss
        merged_df = pd.merge(
            orig_df,
            self.tss_df,
            on=DFSchema.TARGET_GENE,
            how="left",
            suffixes=("", "_ref"),
        )
        orig_df[DFSchema.TARGET_GENE_TSS] = merged_df[DFSchema.START + "_ref"]
        return orig_df


class PredDFLoader(DFLoader):
    def __init__(
        self, filename_tsv: str, tss_filename_bed: str, abc_threshold: float
    ) -> None:
        super().__init__(filename_tsv, tss_filename_bed)
        self.abc_threshold = abc_threshold

    def load(self) -> pd.DataFrame:
        self.df = DFSchema.schematize_pred_df(self.df, self.abc_threshold)
        self.df = self._standardize_gene_tss(self.df)
        return self.df


class CrisprDFLoader(DFLoader):
    def __init__(self, filename_tsv: str, tss_filename_bed: str) -> None:
        super().__init__(filename_tsv, tss_filename_bed)

    def load(self) -> pd.DataFrame:
        self.df = DFSchema.schematize_crispr_df(self.df)
        self.df = self._standardize_gene_tss(self.df)
        return self.df
