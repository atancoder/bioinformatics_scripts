import csv

import pandas as pd
from schema import DFSchema


class DFLoader:
    def __init__(
        self, crispr_filename_tsv: str, pred_filename_tsv: str, tss_filename_bed: str
    ) -> None:
        self.crispr_df = pd.read_csv(
            crispr_filename_tsv, delimiter="\t", compression="gzip"
        )
        self.pred_df = pd.read_csv(
            pred_filename_tsv, delimiter="\t", compression="gzip"
        )
        self.tss_df = self.read_bed_file(tss_filename_bed)

    @classmethod
    def read_bed_file(cls, filename: str) -> pd.DataFrame:
        return pd.read_csv(
            filename,
            sep="\t",
            header=None,
            names=[
                DFSchema.chrom,
                DFSchema.start,
                DFSchema.end,
                DFSchema.target_gene,
            ],
            usecols=range(4),
        )

    def process_dfs(self, abc_threshold: float) -> None:
        self.crispr_df = DFSchema.schematize_crispr_df(self.crispr_df)
        self.pred_df = DFSchema.schematize_pred_df(self.pred_df, abc_threshold)

        self.crispr_df = self.standardize_gene_tss(self.crispr_df)
        self.pred_df = self.standardize_gene_tss(self.pred_df)

    def standardize_gene_tss(self, orig_df: pd.DataFrame) -> pd.DataFrame:
        # filter out rows with genes not in tss_df
        orig_df = orig_df[
            orig_df[DFSchema.target_gene].isin(self.tss_df[DFSchema.target_gene])
        ].reset_index()

        # Use self.tss_df for the target_gene_tss
        merged_df = pd.merge(
            orig_df,
            self.tss_df,
            on=DFSchema.target_gene,
            how="left",
            suffixes=("", "_ref"),
        )
        orig_df[DFSchema.target_gene_tss] = merged_df[DFSchema.start + "_ref"]
        return orig_df
