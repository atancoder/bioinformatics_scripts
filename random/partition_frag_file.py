"""
Partitions a frag file into multiple files based on fragment length

py partition_frag_file.py /oak/stanford/groups/engreitz/Projects/scE2G/data/xu_et_al/fragment.tsv /oak/stanford/groups/engreitz/Users/atan5133/data/scATAC/bins/
"""

import os

import click
import pandas as pd

BINS = [(0, 147), (180, 314), (315, 473), (474, 615)]


@click.command()
@click.argument("frag_file")
@click.argument("out_dir")
def main(frag_file, out_dir):
    frag_basename = os.path.basename(frag_file).split(".")[0]
    frag_df = pd.read_csv(
        frag_file,
        sep="\t",
        names=["chrom", "start", "end"],
        usecols=range(3),
    )
    frag_df["length"] = frag_df["end"] - frag_df["start"]
    for i in range(len(BINS)):
        bin = BINS[i]
        matching_frags = frag_df[
            (frag_df["length"] >= bin[0]) & (frag_df["length"] <= bin[1])
        ]
        bin_file = os.path.join(out_dir, f"{frag_basename}_bin{i}.bed")
        matching_frags.to_csv(bin_file, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
