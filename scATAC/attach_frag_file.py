import click
import pandas as pd
import glob
import os

def find_fragment_file(folder, tissue_name):
    matching_files = glob.glob(os.path.join(folder, f"*{tissue_name}*"))
    if len(matching_files) > 1:
        raise Exception(f"Found {len(matching_files)} matching files with {tissue_name}: {matching_files}")
    elif not matching_files:
        print(f"Couldn't find matching file for tissue: {tissue_name}")
        return None
    return matching_files[0]

def fill_in_fragment_file(metadata_df, fragment_folder):
    known_mappings = {}
    for idx, row in metadata_df.iterrows():
        tissue = row["tissue"]
        if tissue not in known_mappings:
            known_mappings[tissue] = find_fragment_file(fragment_folder, tissue)
        metadata_df.loc[idx, "frag_file"] = known_mappings[tissue]
            
def print_missing_files(frag_files, fragment_folder):
    matching_files = glob.glob(os.path.join(fragment_folder,'*'))
    missing = [m for m in matching_files if m not in set(frag_files)]
    print(f"{missing} files not found in cell metadata")

@click.command
@click.option("--fragments", "fragment_folder", required=True)
@click.option("--metadata", "metadata_tsv", required=True)
@click.option("--output", required=True)
def main(fragment_folder, metadata_tsv, output) -> None:
    metadata_df = pd.read_csv(metadata_tsv, sep="\t")
    fill_in_fragment_file(metadata_df, fragment_folder)
    frag_files = list(metadata_df[~metadata_df["frag_file"].isna()]["frag_file"].drop_duplicates())
    print_missing_files(frag_files, fragment_folder)
    metadata_df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main()