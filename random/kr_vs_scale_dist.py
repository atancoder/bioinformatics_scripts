import hicstraw
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

OUTPUT_FILE = "kr_vs_scale_dist.pdf"
NORM_TYPE = {
    "KR_encode": "https://www.encodeproject.org/files/ENCFF718AWL/@@download/ENCFF718AWL.hic",
    "KR_ncbi": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1551nnn/GSM1551620/suppl/GSM1551620%5FHIC071%5F30.hic",
    "SCALE_K562": "https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic",
    "SCALE_HCT116": "https://www.encodeproject.org/files/ENCFF317OIA/@@download/ENCFF317OIA.hic",
    "SCALE_B_Cell": "https://www.encodeproject.org/files/ENCFF076LWH/@@download/ENCFF076LWH.hic",
    "SCALE_megamap": "https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter2/tissues.hic",
}


def plot_distribution(arr, title, x_label):
    plt.clf()  # make sure we're starting with a fresh plot
    ax = sns.histplot(arr, kde=True, bins=50, kde_kws=dict(cut=3))
    ax.set_title(title)
    ax.set_xlabel(x_label)

    mean = round(np.mean(arr), 3)  # 3 decimal places
    median = round(np.median(arr), 3)
    plt.axvline(x=mean, color="red", linestyle="-", label=f"Mean={mean}")
    plt.axvline(x=median, color="green", linestyle="-", label=f"Median={median}")
    plt.legend()

    return ax.get_figure()


with PdfPages(OUTPUT_FILE) as pdf_writer:
    for norm in NORM_TYPE:
        if "KR" in norm:
            chr = "22"
            norm_type = "KR"
        else:
            chr = "chr22"
            norm_type = "SCALE"
        records = hicstraw.straw(
            "norm",
            norm_type,
            NORM_TYPE[norm],
            chr,
            chr,
            "BP",
            5000,
        )
        counts = np.array([r.counts for r in records])
        below_threshold = len(counts[counts<.25])
        pct = round(below_threshold / len(counts), 5) * 100
        print(f"{norm} has {below_threshold}/{len(counts)} below the .25 threshold. This is {pct}%")
        norm_counts = counts[counts<3]
        figure = plot_distribution(
            norm_counts, f"{norm} norms for chrom: {chr}", "Value"
        )
        pdf_writer.savefig(figure)
