import csv
import os
import time
from typing import Dict, List, NamedTuple, Set, Tuple

import pandas as pd
from df_loader import DFLoader
from overlaps import (
    compute_all_overlaps,
    read_overlaps_from_file,
    write_overlaps_to_file,
)
from schema import DFSchema

CRISPR_FILENAME = "resources/example/EPCrisprBenchmark_Fulco2019_K562_GRCh38.tsv.gz"
PRED_FILENAME = "resources/example/ABC_K562_Fulco2019Genes_GRCh38.tsv.gz"
OVERLAP_FILENAME = "overlaps.csv"
ABC_THRESHOLD = 0.02
TSS_REF_FILE = "resources/genome_annotations/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed"


def main():
    if os.path.exists(OVERLAP_FILENAME):
        overlaps = read_overlaps_from_file(OVERLAP_FILENAME)
    else:
        start = time.time()
        df_loader = DFLoader(CRISPR_FILENAME, PRED_FILENAME, TSS_REF_FILE)
        df_loader.process_dfs(ABC_THRESHOLD)
        crispr_df, pred_df = df_loader.crispr_df, df_loader.pred_df
        overlaps = compute_all_overlaps(crispr_df, pred_df)
        write_overlaps_to_file(overlaps, OVERLAP_FILENAME)
        print(f"Found and saved overlaps in {time.time() - start} seconds")


if __name__ == "__main__":
    main()
