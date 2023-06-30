from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from schema import DFSchema


class CategoryLabeler(ABC):
    COL_SUFFIX = "_Category"

    @classmethod
    def name(cls):
        # e.g FalsePos_Category
        if cls is CategoryLabeler:
            raise Exception("Calling name() on Abstract Class")
        return cls.__name__ + cls.COL_SUFFIX

    @classmethod
    @abstractmethod
    def get_values(cls, df: pd.DataFrame) -> pd.Series:
        """
        Values of the category column
        e.g pd.Series(["false_pos", "false_pos", "false_neg", ...])
        """
        raise NotImplementedError()

    @classmethod
    def create_category(cls, df: pd.DataFrame) -> None:
        df[cls.name()] = cls.get_values(df)

    @classmethod
    def summarize_category_count(cls, df: pd.DataFrame) -> pd.Series:
        column_name = cls.name()
        col = df[column_name]
        return col.value_counts()


class FalsePos(CategoryLabeler):
    @classmethod
    def get_values(cls, df: pd.DataFrame) -> pd.Series:
        return (df[DFSchema.is_significant + DFSchema.crispr_suffix] == False) & (
            df[DFSchema.is_significant + DFSchema.pred_suffix] == True
        )

    @classmethod
    def summarize_category_count(cls, df: pd.DataFrame) -> str:
        return super().summarize_category_count(df).filter([True])


class FalseNeg(CategoryLabeler):
    @classmethod
    def get_values(cls, df: pd.DataFrame) -> pd.Series:
        return (df[DFSchema.is_significant + DFSchema.crispr_suffix] == True) & (
            df[DFSchema.is_significant + DFSchema.pred_suffix] == False
        )

    @classmethod
    def summarize_category_count(cls, df: pd.DataFrame) -> str:
        return super().summarize_category_count(df).filter([True])


class DistToTSSCategory(CategoryLabeler):
    SMALL = 1e5  # 100KB
    MEDIUM = 1e6  # 1 MB

    @classmethod
    def categorize_distance(cls, abs_dist_to_tss: int) -> str:
        if abs_dist_to_tss <= cls.SMALL:
            return "small"
        elif abs_dist_to_tss <= cls.MEDIUM:
            return "medium"
        else:
            return "large"

    @classmethod
    def get_values(cls, df: pd.DataFrame) -> pd.Series:
        dist_to_tss = np.abs(df[DFSchema.dist_to_tss + DFSchema.crispr_suffix])
        return dist_to_tss.apply(cls.categorize_distance)
