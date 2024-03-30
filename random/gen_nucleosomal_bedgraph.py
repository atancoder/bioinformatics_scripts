import click
import pandas as pd

"""
py gen_nucleosomal_bedgraph.py /oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_ATAC/binned_atac/ENCFF534DCE_bin0_abc.bedgraph /oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_ATAC/binned_atac/ENCFF534DCE_bin1_abc.bedgraph /oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_ATAC/binned_atac/ENCFF534DCE_bin2_abc.bedgraph /oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_ATAC/binned_atac/ENCFF534DCE_bin3_abc.bedgraph /oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_ATAC/binned_atac/ENCFF534DCE_NR_abc.bedgraph
"""


@click.command
@click.argument("bin0")
@click.argument("bin1")
@click.argument("bin2")
@click.argument("bin3")
@click.argument("output")
def main(bin0, bin1, bin2, bin3, output):
    bin0_bed = pd.read_csv(
        bin0,
        sep="\t",
        names=["chrom", "start", "end", "count"],
    )
    bin1_bed = pd.read_csv(
        bin1,
        sep="\t",
        names=["chrom", "start", "end", "count"],
    )
    bin2_bed = pd.read_csv(
        bin2,
        sep="\t",
        names=["chrom", "start", "end", "count"],
    )
    bin3_bed = pd.read_csv(
        bin3,
        sep="\t",
        names=["chrom", "start", "end", "count"],
    )

    # Apply multipliers
    bin0_bed["count"] *= -1
    bin1_bed["count"] *= 1
    bin2_bed["count"] *= 2
    bin3_bed["count"] *= 3

    nucleosomal_bedgraph = bin0_bed[["chrom", "start", "end"]]
    nucleosomal_bedgraph["count"] = (
        bin0_bed["count"] + bin1_bed["count"] + bin2_bed["count"] + bin3_bed["count"]
    )
    nucleosomal_bedgraph["count"] = nucleosomal_bedgraph["count"].apply(
        lambda x: max(x, 0)
    )
    nucleosomal_bedgraph.to_csv(output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
