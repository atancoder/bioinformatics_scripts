import glob
import os
import re

import click
import pandas as pd

"""
Fill in the metadata df with new column that references the fragment file
"""

def find_fragment_file(folder, sample):
    matching_files = glob.glob(os.path.join(folder, f"*{sample}*"))
    if len(matching_files) > 1:
        raise Exception(f"Found {len(matching_files)} matching files with {sample}: {matching_files}")
    elif not matching_files:
        print(f"Couldn't find matching file for sample: {sample}")
        return None
    return matching_files[0]

def fill_in_fragment_file(metadata_df, fragment_folder):
    known_mappings = {}
    for idx, row in metadata_df.iterrows():
        sample = row["sample"]
        if sample not in known_mappings:
            known_mappings[sample] = find_fragment_file(fragment_folder, sample)
        metadata_df.loc[idx, "frag_file"] = known_mappings[sample]
            
def print_missing_files(frag_files, fragment_folder):
    matching_files = glob.glob(os.path.join(fragment_folder,'*'))
    missing = [m for m in matching_files if m not in set(frag_files)]
    print(f"{missing} files not found in cell metadata")

def sanitize_filename(filename):
    """
    We don't want spaces or parenthesis for filenames
    """
    return re.sub(r'[^\w]', '_', filename)


def sanitize_cell_type(df):
    num_cell_types = len(df["cell_type"].drop_duplicates())
    new_names = df["cell_type"].apply(sanitize_filename)
    if num_cell_types != len(new_names.drop_duplicates()):
        diff = len(new_names.drop_duplicates()) - num_cell_types
        raise Exception(f"Sanitizing led to a reduction of {diff} cell types")
    df["cell_type"] = new_names    

def rename_columns(df):
    # Standardize column name
    df.rename(columns={"cell type": "cell_type"}, inplace=True)
    df.rename(columns={"tissue": "sample"}, inplace=True)

@click.command
@click.option("--fragments", "fragment_folder", required=True)
@click.option("--metadata", "metadata_tsv", required=True)
@click.option("--output", required=True)
def main(fragment_folder, metadata_tsv, output) -> None:
    metadata_df = pd.read_csv(metadata_tsv, sep="\t")
    rename_columns(metadata_df)
    sanitize_cell_type(metadata_df)
    fill_in_fragment_file(metadata_df, fragment_folder)
    frag_files = list(metadata_df[~metadata_df["frag_file"].isna()]["frag_file"].drop_duplicates())
    print_missing_files(frag_files, fragment_folder)
    metadata_df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main()