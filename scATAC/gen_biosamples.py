import os
import glob
import click
import pandas as pd

CONFIG_TEMPLATE_FILE = "https://raw.githubusercontent.com/broadinstitute/ABC-Enhancer-Gene-Prediction/dev/config/config_biosamples_template.tsv"


@click.command()
@click.option("--tagAlign_folder", type=str, required=True)
@click.option("--output", type=str, required=True)
def main(tagalign_folder, output):
    biosamples = pd.read_csv(CONFIG_TEMPLATE_FILE, sep="\t")
    tagAlign_files = glob.glob(os.path.join(tagalign_folder, "*tagAlign.gz"))
    for i, tagAlign in enumerate(tagAlign_files):
        base_file = os.path.basename(tagAlign)
        cell_type = base_file[: base_file.index(".tagAlign.gz")]
        biosample = {
            "biosample": cell_type,
            "ATAC": tagAlign,
            "default_accessibility_feature": "ATAC",
        }
        biosamples.loc[i] = biosample

    biosamples.to_csv(output, sep="\t", index=False)
    print(f"Generated config: {output}")


if __name__ == "__main__":
    main()
