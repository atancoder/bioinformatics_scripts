import csv
import os
import time
from typing import Dict, List, NamedTuple, Set, Tuple

import bioframe as bf
import pandas as pd
from schema import DFSchema


def compute_all_overlaps(
    crispr_df: pd.DataFrame, pred_df: pd.DataFrame
) -> pd.DataFrame:
    cols = [DFSchema.chrom, DFSchema.start, DFSchema.end]
    overlap_df = bf.overlap(
        crispr_df,
        pred_df,
        how="inner",
        return_index=True,
        cols1=cols,
        cols2=cols,
        on=[DFSchema.target_gene],
        suffixes=(DFSchema.crispr_suffix, DFSchema.pred_suffix),
    )
    return overlap_df


def write_overlaps_to_file(overlap_df: pd.DataFrame, file: str) -> None:
    overlap_df.to_csv(file, index=False)


def read_overlaps_from_file(file: str) -> pd.DataFrame:
    return pd.read_csv(file)
