import click
import metaseq
import numpy as np
import pybedtools

REGION_SLOP = 2000
NUM_BINS = 400


def get_tss_enrichment(
    input_file,
    tss,
    chrom_sizes,
    bins=NUM_BINS,
    bp_edge=REGION_SLOP,
    processes=8,
):
    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chrom_sizes)

    # Load the bam file
    # Need to shift reads and just get ends, just load bed file?
    if input_file.endswith("bam"):
        signal = metaseq.genomic_signal(input_file, "bam")
    else:
        signal = metaseq.genomic_signal(input_file, "bed")

    signal_array = signal.array(
        tss_ext,
        bins=bins,
        # processes=processes,
        stranded=True,
    )
    # Normalization (Greenleaf style): Find the avg height
    # at the end bins and take fold change over that
    # Use enough bins to cover 100 bp on either end
    num_edge_bins = int(100 / (2 * bp_edge / bins))
    edge_signal = np.concatenate((signal_array[:, :num_edge_bins], signal_array[:, -num_edge_bins:]), axis=1)
    avg_noise = edge_signal.sum(axis=1) / (2 * num_edge_bins)
    min_noise = avg_noise[avg_noise!=0].min()
    avg_noise[avg_noise==0] = min_noise
    signal_array /= avg_noise.reshape(-1, 1)

    highest_tss_scores = signal_array.max(axis=1)
    return highest_tss_scores.mean()


@click.command()
@click.option("-i", "inputs", help="Can be bam file or fragment file. Can be comma separated")
@click.option("--tss", help="bed file consisting of TSS cordinates of 1bp")
@click.option("--chrom_sizes")
def main(inputs, tss, chrom_sizes):
    for input in inputs.split(","):
        input = input.strip()
        tss_score = get_tss_enrichment(input, tss, chrom_sizes)
        print("TSS Score for {}: {}".format(input, tss_score))


"""
conda activate encd-atac-py2
py tss_enrichment.py -i /oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_DNASE/ENCFF860XAE.sorted.se.bam,/oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_DNASE/duke/ENCFF987IUK_sorted.bam,/oak/stanford/groups/engreitz/Users/atan5133/data/scATAC/xu_150bp.tagAlign.gz,/oak/stanford/groups/engreitz/Projects/scE2G/data/xu_et_al/xu_K562_sorted.tagAlign.gz,/oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_ATAC/ENCFF512VEZ.sorted.bam,/oak/stanford/groups/engreitz/Users/atan5133/data/ENCODE/K562_ATAC/ENCFF534DCE.sorted.bam --tss chr22_tss.bed --chrom_sizes /oak/stanford/groups/engreitz/Users/atan5133/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.chrom.sizes.tsv
"""

if __name__ == "__main__":
    main()
