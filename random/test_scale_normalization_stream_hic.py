import os

import hicstraw
import numpy as np
import pandas as pd
import scipy.sparse as ssp

HIC_FILE = "https://www.encodeproject.org/files/ENCFF309UNV/@@download/ENCFF309UNV.hic"
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


def get_chrom_format(hic: hicstraw.HiCFile, chromosome):
    """
    hic files can have 'chr1' or just '1' as the chromosome name
    we need to make sure we're using the format consistent with
    the hic file
    """
    hic_chrom_names = [chrom.name for chrom in hic.getChromosomes()]
    if hic_chrom_names[1].startswith("chr"):  # assume index 1 should be chr1
        return chromosome
    else:
        return chromosome[3:]


def get_hic_mat(hic_file, chromosome, resolution=5000):
    hic = hicstraw.HiCFile(hic_file)
    chromosome = get_chrom_format(hic, chromosome)
    results = hicstraw.straw(
        "observed",
        "SCALE",
        hic_file,
        f"{chromosome}:0:1000000",
        f"{chromosome}:0:1000000",
        "BP",
        resolution,
    )
    results = results[:10]
    for result in results:
        print(result.binX, result.binY, result.counts)


hic_mat = get_hic_mat(HIC_FILE, "chr1")
# Check similarity to (base) [atan5133@sh02-ln03 login /oak/stanford/groups/engreitz/Projects/ABC/HiC/w80-mucosa-of-descending-colon-intact/chr1]$ zcat chr1.INTERSCALEobserved.gz | head
