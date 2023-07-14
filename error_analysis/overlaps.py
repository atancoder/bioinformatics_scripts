import csv
import os
import time
from typing import Dict, List, NamedTuple, Set, Tuple

import bioframe as bf
import pandas as pd

from schema import DFSchema

OVERLAP_SUFFIXES = ["_a", "_b"]


def compute_all_overlaps(
    df_a: pd.DataFrame, df_b: pd.DataFrame, suffixes=OVERLAP_SUFFIXES
) -> pd.DataFrame:
    cols = [DFSchema.CHROM, DFSchema.START, DFSchema.END]
    overlap_df = bf.overlap(
        df_a,
        df_b,
        how="inner",
        return_index=True,
        cols1=cols,
        cols2=cols,
        on=[DFSchema.TARGET_GENE],
        suffixes=suffixes,
    )
    return overlap_df


def write_overlaps_to_file(overlap_df: pd.DataFrame, file: str) -> None:
    overlap_df.to_csv(file, index=False)


def read_overlaps_from_file(file: str) -> pd.DataFrame:
    return pd.read_csv(file)
