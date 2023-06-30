from typing import Dict

import pandas as pd


class DFSchema:
    """
    Both Dataframes will adopt this schema
    """

    chrom = "chrom"
    start = "start"
    end = "end"
    target_gene = "TargetGene"
    target_gene_tss = "TargetGeneTSS"
    is_significant = "IsSignficant"
    dist_to_tss = "DistanceToTSS"

    crispr_suffix = "_crispr"
    pred_suffix = "_pred"

    @classmethod
    def get_crispr_schema_map(cls) -> Dict[str, str]:
        return {
            "chrom": cls.chrom,
            "chromStart": cls.start,
            "chromEnd": cls.end,
            "measuredGeneSymbol": cls.target_gene,
            "startTSS": cls.target_gene_tss,
            "Significant": cls.is_significant,
        }

    @classmethod
    def get_pred_schema_map(cls) -> Dict[str, str]:
        return {
            "chr": cls.chrom,
            "start": cls.start,
            "end": cls.end,
            "TargetGene": cls.target_gene,
            "TargetGeneTSS": cls.target_gene_tss,
        }

    @classmethod
    def schematize_crispr_df(cls, df: pd.DataFrame) -> pd.DataFrame:
        schema_map = cls.get_crispr_schema_map()
        df = df.rename(columns=schema_map)
        cls.add_dist_to_gene_col(df)
        return df

    @classmethod
    def schematize_pred_df(cls, df: pd.DataFrame, threshold: float) -> pd.DataFrame:
        schema_map = cls.get_pred_schema_map()
        df = df.rename(columns=schema_map)
        cls.add_dist_to_gene_col(df)
        df[cls.is_significant] = df["powerlaw.Score"] >= threshold
        return df

    @classmethod
    def add_dist_to_gene_col(cls, df: pd.DataFrame) -> None:
        df[cls.dist_to_tss] = df[cls.target_gene_tss] - df[cls.start]
