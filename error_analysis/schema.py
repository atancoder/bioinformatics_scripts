from typing import Dict

import pandas as pd


class DFSchema:
    """
    Both Dataframes will adopt this schema
    """

    CHROM = "chrom"
    START = "start"
    END = "end"
    TARGET_GENE = "TargetGene"
    TARGET_GENE_TSS = "TargetGeneTSS"
    IS_SIGNIFICANT = "IsSignificant"
    DIST_TO_TSS = "DistanceToTSS"

    CRISPR_SUFFIX = "_crispr"
    PRED_SUFFIX = "_pred"

    FROM_MULTIPLE_PREDICTIONS = "from_mult_pred"

    @classmethod
    def get_crispr_schema_map(cls) -> Dict[str, str]:
        return {
            "chrom": cls.CHROM,
            "chromStart": cls.START,
            "chromEnd": cls.END,
            "measuredGeneSymbol": cls.TARGET_GENE,
            "startTSS": cls.TARGET_GENE_TSS,
            "Significant": cls.IS_SIGNIFICANT,
        }

    @classmethod
    def get_pred_schema_map(cls) -> Dict[str, str]:
        return {
            "chr": cls.CHROM,
            "start": cls.START,
            "end": cls.END,
            "TargetGene": cls.TARGET_GENE,
            "TargetGeneTSS": cls.TARGET_GENE_TSS,
        }

    @classmethod
    def schematize_crispr_df(cls, df: pd.DataFrame) -> pd.DataFrame:
        schema_map = cls.get_crispr_schema_map()
        df = df.rename(columns=schema_map)
        cls.add_dist_to_gene_col(df)
        return df

    @classmethod
    def schematize_pred_df(
        cls, df: pd.DataFrame, threshold: float, score_col: str
    ) -> pd.DataFrame:
        schema_map = cls.get_pred_schema_map()
        df = df.rename(columns=schema_map)
        cls.add_dist_to_gene_col(df)
        df[cls.IS_SIGNIFICANT] = df[score_col] >= threshold
        return df

    @classmethod
    def add_dist_to_gene_col(cls, df: pd.DataFrame) -> None:
        df[cls.DIST_TO_TSS] = df[cls.TARGET_GENE_TSS] - df[cls.START]
