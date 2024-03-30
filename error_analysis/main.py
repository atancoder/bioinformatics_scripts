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
PRED_FILENAME = "/oak/stanford/groups/engreitz/Users/atan5133/abc_run_comparisons/current_ABC_results/K562_dhs_pl/Predictions/EnhancerPredictionsAllPutative.tsv.gz"
OVERLAP_FILENAME = "crispr_atac_pred_overlaps.tsv"
ABC_THRESHOLD = 0.016
TSS_REF_FILE = "resources/genome_annotations/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed"


def main():
    if os.path.exists(OVERLAP_FILENAME):
        overlaps = read_overlaps_from_file(OVERLAP_FILENAME)
    else:
        start = time.time()
        pred_df = PredDFLoader(PRED_FILENAME, TSS_REF_FILE, ABC_THRESHOLD).load()
        crispr_df = CrisprDFLoader(CRISPR_FILENAME, TSS_REF_FILE).load()
        overlaps = compute_crispr_overlaps(crispr_df, pred_df, ABC_THRESHOLD)
        write_overlaps_to_file(overlaps, OVERLAP_FILENAME)
        print(f"Found and saved overlaps in {time.time() - start} seconds")


if __name__ == "__main__":
    main()
