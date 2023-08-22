import csv
import os
import time
from typing import Dict, List, NamedTuple, Set, Tuple
import numpy as np

import bioframe as bf
import pandas as pd

from schema import DFSchema

OVERLAP_SUFFIXES = ["_a", "_b"]


def compute_all_overlaps(
    df_a: pd.DataFrame, df_b: pd.DataFrame, suffixes=OVERLAP_SUFFIXES
) -> pd.DataFrame:
    """
    Compure overlaps between 2 DFs via inner join
    """
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

def compute_crispr_overlaps(
    crispr_df: pd.DataFrame, pred_df: pd.DataFrame, 
) -> pd.DataFrame:
    """
    For each E-G pair in CRISPR df, finds the entry from pred_df that overlaps it
    There may not be an overlapping prediction

    If there are multiple pred matches to a CRISPR DF, we will aggregate those matches and present
    a score (ABC, activity, contact) in a later step (merge_multiple_predictions)
    """
    cols = [DFSchema.CHROM, DFSchema.START, DFSchema.END]
    overlap_df = bf.overlap(
        crispr_df,
        pred_df,
        how="left",
        cols1=cols,
        cols2=cols,
        on=[DFSchema.TARGET_GENE],
        suffixes=[DFSchema.CRISPR_SUFFIX, DFSchema.PRED_SUFFIX],
    )
    overlap_df[DFSchema.IS_SIGNIFICANT + DFSchema.PRED_SUFFIX].fillna(False, inplace=True)
    return overlap_df

def merge_multiple_predictions(overlap_df: pd.DataFrame, agg_fn=np.mean) -> pd.DataFrame:
    """
    When there are multiple matching predictions for a crispri experiment, 
    we aggregate the scores (default mean)
    """
    name_col = "name" + DFSchema.CRISPR_SUFFIX
    duplicated_names = overlap_df[overlap_df.duplicated(name_col)][name_col]
    new_df = overlap_df.drop_duplicates(subset=[name_col]).reset_index(drop=True)
    new_df[DFSchema.FROM_MULTIPLE_PREDICTIONS] = False

    for name in duplicated_names:
        entries = overlap_df[overlap_df[name_col] == name]
        new_score = agg_fn(entries["ABC.Score" + DFSchema.PRED_SUFFIX])
        new_entry = new_df[new_df[name_col] == name].index  # Should only be 1 b/c we removed duplicates
        new_df.loc[new_entry, "ABC.Score" + DFSchema.PRED_SUFFIX] = new_score
        new_df.loc[new_entry, DFSchema.FROM_MULTIPLE_PREDICTIONS] = True
    return new_df


def write_overlaps_to_file(overlap_df: pd.DataFrame, file: str) -> None:
    overlap_df.to_csv(file, index=False)


def read_overlaps_from_file(file: str) -> pd.DataFrame:
    return pd.read_csv(file)
