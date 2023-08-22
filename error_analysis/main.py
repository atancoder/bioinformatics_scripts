import csv
import os
import time
from typing import Dict, List, NamedTuple, Set, Tuple

import pandas as pd
from df_loader import CrisprDFLoader, PredDFLoader
from overlaps import (
    compute_crispr_overlaps,
    read_overlaps_from_file,
    write_overlaps_to_file,
)
from schema import DFSchema

CRISPR_FILENAME = "/oak/stanford/groups/engreitz/Projects/Benchmarking/CRISPR_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"
PRED_FILENAME = "/oak/stanford/groups/engreitz/Users/rosaxma/ABC_single_cell_hg38/K562_ABC_pipeline_benchmark/ENCODE_ABC_output/Weilin_filtered_frag/Predictions/EnhancerPredictionsAllPutative.txt.gz"
OVERLAP_FILENAME = "crispr_pred_overlaps.csv"
ABC_THRESHOLD = 0.02
TSS_REF_FILE = "resources/genome_annotations/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed"


def main():
    if os.path.exists(OVERLAP_FILENAME):
        overlaps = read_overlaps_from_file(OVERLAP_FILENAME)
    else:
        start = time.time()
        pred_df = PredDFLoader(PRED_FILENAME, TSS_REF_FILE, ABC_THRESHOLD).load()
        crispr_df = CrisprDFLoader(CRISPR_FILENAME, TSS_REF_FILE).load()
        overlaps = compute_crispr_overlaps(
            crispr_df, pred_df
        )
        write_overlaps_to_file(overlaps, OVERLAP_FILENAME)
        print(f"Found and saved overlaps in {time.time() - start} seconds")


if __name__ == "__main__":
    main()
