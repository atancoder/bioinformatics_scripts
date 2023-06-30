import inspect
from abc import ABC, abstractmethod
from typing import Dict, List, NamedTuple, Set, Tuple, Type

import pandas as pd
from category_labelers import CategoryLabeler, DistToTSSCategory, FalseNeg, FalsePos
from schema import DFSchema

LABELERS: List[Type[CategoryLabeler]] = [FalsePos, FalseNeg, DistToTSSCategory]


def get_misclassifications(overlap_df: pd.DataFrame) -> pd.DataFrame:
    return overlap_df[
        overlap_df[DFSchema.is_significant + DFSchema.crispr_suffix]
        != overlap_df[DFSchema.is_significant + DFSchema.pred_suffix]
    ].reset_index()


def label_misclassifications(df: pd.DataFrame) -> None:
    for labeler in LABELERS:
        labeler.create_category(df)
