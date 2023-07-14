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
        return (df[DFSchema.IS_SIGNIFICANT + DFSchema.CRISPR_SUFFIX] == False) & (
            df[DFSchema.IS_SIGNIFICANT + DFSchema.PRED_SUFFIX] == True
        )

    @classmethod
    def summarize_category_count(cls, df: pd.DataFrame) -> str:
        return super().summarize_category_count(df).filter([True])


class FalseNeg(CategoryLabeler):
    @classmethod
    def get_values(cls, df: pd.DataFrame) -> pd.Series:
        return (df[DFSchema.IS_SIGNIFICANT + DFSchema.CRISPR_SUFFIX] == True) & (
            df[DFSchema.IS_SIGNIFICANT + DFSchema.PRED_SUFFIX] == False
        )

    @classmethod
    def summarize_category_count(cls, df: pd.DataFrame) -> str:
        return super().summarize_category_count(df).filter([True])


class DistToTSSSize(CategoryLabeler):
    SMALL = 1e4  # 10KB
    MEDIUM = 1e5  # 100KB

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
        dist_to_tss = np.abs(df[DFSchema.DIST_TO_TSS + DFSchema.CRISPR_SUFFIX])
        return dist_to_tss.apply(cls.categorize_distance)


class Top5Gene(CategoryLabeler):
    # Class works a bit different than other labelers because
    # we don't need to create a new column

    @classmethod
    def get_values(cls, df: pd.DataFrame) -> pd.Series:
        """
        Values of the category column
        e.g pd.Series(["false_pos", "false_pos", "false_neg", ...])
        """
        raise NotImplementedError()

    @classmethod
    def create_category(cls, df: pd.DataFrame) -> None:
        return

    @classmethod
    def summarize_category_count(cls, df: pd.DataFrame) -> pd.Series:
        name = cls.name()
        top_5 = df[DFSchema.TARGET_GENE + DFSchema.CRISPR_SUFFIX].value_counts()[:5]
        top_5 = top_5.rename_axis(name)
        return top_5
