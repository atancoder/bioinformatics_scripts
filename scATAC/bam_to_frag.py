import os

import click
import pysam


def bam_to_frag(in_path, out_path, shift_plus, shift_minus, is_sc=True):
    """
    Convert coordinate-sorted BAM file to a fragment file format, while adding Tn5 coordinate adjustment
    BAM should be pre-filtered for PCR duplicates, secondary alignments, and unpaired reads
    Output fragment file is sorted by chr, start, end, barcode
    """

    input = pysam.AlignmentFile(in_path, "rb")
    with open(out_path, "w") as out_file:
        buf = []
        curr_pos = None
        for read in input:
            if read.flag & 16 == 16:
                continue  # ignore reverse (coordinate-wise second) read in pair

            chromosome = read.reference_name
            tlen = read.template_length
            if tlen <= 0:
                continue

            start = read.reference_start + shift_plus
            end = read.reference_start + tlen + shift_minus
            if is_sc:
                cell_barcode = read.get_tag("CB")
                data = (chromosome, start, end, cell_barcode, 1)
            else:
                data = (chromosome, start, end)
            pos = (chromosome, start)

            if pos == curr_pos:
                buf.append(data)
            else:
                buf.sort()
                for i in buf:
                    print(*i, sep="\t", file=out_file)
                buf.clear()
                buf.append(data)
                curr_pos = pos


@click.command()
@click.option("--bam")
def main(bam):
    filename = os.path.basename(bam).split(".bam")[0]
    frag_file = os.path.join(os.path.dirname(bam), f"{filename}.tsv")
    bam_to_frag(bam, frag_file, 4, -5, is_sc=False)


if __name__ == "__main__":
    main()
