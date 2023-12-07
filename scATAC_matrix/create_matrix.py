import bioframe as bf
import click
import pandas as pd
import time
from collections import defaultdict
import warnings
import matplotlib.pyplot as plt

CRISPR_FILENAME = "/oak/stanford/groups/engreitz/Projects/Benchmarking/CRISPR_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"
TSS_FILE = "/oak/stanford/groups/engreitz/Users/atan5133/abc_run_comparisons/ABC-Enhancer-Gene-Prediction/reference/hg38/CollapsedGeneBounds.hg38.TS500bp.bed"
FRAGMENT_FILE = "/oak/stanford/groups/engreitz/Users/atan5133/scATAC_data/fragments_rosa_chr8.tsv"
GENE = "MYC"

def get_enh_regions(crispr_df, chrom):
    df = crispr_df[crispr_df["chrom"] == chrom]
    return crispr_df[["chrom", "start", "end"]].drop_duplicates().reset_index(drop=True)

def create_accessibility_matrix(overlap_df):
    barcodes = overlap_df["barcode_"].unique()
    max_enh_idx = overlap_df["index"].max()
    matrix = pd.DataFrame(0,index=barcodes, columns=range(max_enh_idx+1))
    for barcode, df in overlap_df.groupby("barcode_"):
        enh_indices = df["index"]
        matrix.loc[barcode] = {idx: 1 for idx in enh_indices}
    matrix = matrix.fillna(0).astype(int)
    return matrix

def main():
    crispr_df = pd.read_csv(CRISPR_FILENAME, sep="\t")
    crispr_df = crispr_df.rename(columns={"chromStart": "start", "chromEnd": "end"})
    tss_df = pd.read_csv(TSS_FILE, sep="\t")
    tss_df = tss_df.rename(columns={"#chr": "chrom"})
    atac_df = pd.read_csv(
        FRAGMENT_FILE,
        sep="\t",
        names=["chrom", "start", "end", "barcode", "readSupport"],
    )
    for chrom in ["chr8"]:
        chr_atac_df = atac_df[atac_df["chrom"] == chrom]
        crispr_enh_df = get_enh_regions(crispr_df, chrom)
        overlap_df = bf.overlap(crispr_enh_df, chr_atac_df, how="inner", return_index=True)
        matrix = create_accessibility_matrix(overlap_df)

        

if __name__ == "__main__":
    main()
