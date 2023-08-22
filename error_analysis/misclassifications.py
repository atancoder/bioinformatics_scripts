import inspect
from abc import ABC, abstractmethod
from typing import Dict, List, NamedTuple, Set, Tuple, Type

import pandas as pd

from category_labelers import *
from schema import DFSchema

LABELERS: List[Type[CategoryLabeler]] = [FalsePos, FalseNeg, DistToTSSSize, Top5Gene, MultiplePredictions]


def get_misclassifications(overlap_df: pd.DataFrame) -> pd.DataFrame:
    return overlap_df[
        overlap_df[DFSchema.IS_SIGNIFICANT + DFSchema.CRISPR_SUFFIX]
        != overlap_df[DFSchema.IS_SIGNIFICANT + DFSchema.PRED_SUFFIX]
    ].reset_index(drop=True)


def label_misclassifications(misclass_df: pd.DataFrame) -> None:
    for labeler in LABELERS:
        labeler.create_category(misclass_df)


def groupby_gene(misclass_df: pd.DataFrame) -> Dict[str, List[int]]:
    return misclass_df.groupby([DFSchema.TARGET_GENE + DFSchema.CRISPR_SUFFIX]).groups
