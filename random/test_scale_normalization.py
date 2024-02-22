import os

import numpy as np
import pandas as pd
import scipy.sparse as ssp

BASE_DIR = "/oak/stanford/groups/engreitz/Projects/ABC/HiC"
DIRS = [
    "K562_intact_Hi-C",
    "Monocytes_intact_Hi-C",
    "W76_mucosa_of_descending_colon",
    "uw040-heart-left-ventricle-intact",
    "w62-left-lung",
    "W80_psoas_muscle_intact_dnase_hic",
]

DIRS = [os.path.join(BASE_DIR, dir) for dir in DIRS]

chromosomes = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
]


def get_hic_mat(obs_name, norm_file, resolution=5000):
    norms = pd.read_csv(norm_file, header=None)
    hic_size = norms.shape[0]

    HiC = pd.read_table(
        obs_name,
        names=["bin1", "bin2", "hic_contact"],
        header=None,
        engine="c",
        memory_map=True,
    )
    row = np.floor(HiC.bin1.values / resolution).astype(int)
    col = np.floor(HiC.bin2.values / resolution).astype(int)
    dat = HiC.hic_contact.values

    mask = row != col  # off-diagonal
    row2 = col[mask]  # note the row/col swap
    col2 = row[mask]
    dat2 = dat[mask]

    # concat and create
    row = np.hstack((row, row2))
    col = np.hstack((col, col2))
    dat = np.hstack((dat, dat2))
    return ssp.csr_matrix((dat, (row, col)), (hic_size, hic_size))

rows = []
for dir in DIRS:
    print(f"Processing {dir}")
    for chr in chromosomes:
        obs_name = os.path.join(dir, chr, f"{chr}.INTERSCALEobserved.gz")
        norm_name = os.path.join(dir, chr, f"{chr}.INTERSCALEnorm.gz")
        hic_mat = get_hic_mat(obs_name, norm_name)
        temp = hic_mat
        temp.data = np.nan_to_num(temp.data, copy=False)
        sums = temp.sum(axis=0)
        sums = sums[~np.isnan(sums)]
        row_mean = np.mean(sums[sums > 0])
        print(f"HiC Matrix {obs_name} has row sums of {row_mean}")
        sums = temp.sum(axis=1)
        sums = sums[~np.isnan(sums)]
        row = {"dir": dir, "chr": chr, "row_mean": row_mean}
        rows.append(row)

df = pd.DataFrame(rows)
df.to_csv("hic_stats.tsv", sep="\t", index=False)
