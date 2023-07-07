import inspect
from abc import ABC, abstractmethod
from typing import Dict, List, NamedTuple, Set, Tuple, Type

import pandas as pd

from category_labelers import *
from schema import DFSchema

LABELERS: List[Type[CategoryLabeler]] = [FalsePos, FalseNeg, DistToTSSSize, Top5Gene]


def get_misclassifications(overlap_df: pd.DataFrame) -> pd.DataFrame:
    return overlap_df[
        overlap_df[DFSchema.is_significant + DFSchema.crispr_suffix]
        != overlap_df[DFSchema.is_significant + DFSchema.pred_suffix]
    ].reset_index()


def label_misclassifications(misclass_df: pd.DataFrame) -> None:
    for labeler in LABELERS:
        labeler.create_category(misclass_df)

def groupby_gene(misclass_df: pd.DataFrame) -> Dict[str, List[int]]:
    return misclass_df.groupby([DFSchema.target_gene + DFSchema.crispr_suffix]).groups