import click
import pandas as pd
import glob
import time
import os
import re
import gzip
import concurrent.futures
import subprocess
from collections import defaultdict

# Num rows to read from csv at a time. We have to be careful
# with how much we read to prevent OOMs
CHUNK_SIZE = 30  

def map_columns_to_cell_type(columns, metadata_df):
    cell_type_cols = {}
    for col in columns:
        first, second, third, fourth = col.split(".")
        cell_id = f"{first}-{second}#{third}-{fourth}"
        cell_type = metadata_df.loc[cell_id]["cell_type"]
        if cell_type not in cell_type_cols:
            cell_type_cols[cell_type] = ["Gene"]
        cell_type_cols[cell_type].append(col)
    
    return cell_type_cols

def create_cell_type_rna_files(rna_df, cell_type_cols, output_folder):
    cell_type_files = {}
    for cell_type, cols in cell_type_cols.items():
        rna_file = os.path.join(output_folder, f"{cell_type}_rna_count_matrix.csv.gz")
        rna_df[cols].to_csv(rna_file, index=False, mode='x')
        cell_type_files[cell_type] = rna_file
    return cell_type_files

@click.command
@click.option("--rna", "rna_tsv", required=True)
@click.option("--metadata", "metadata_tsv", required=True)
@click.option("--output_folder", required=True)
def main(rna_tsv, metadata_tsv, output_folder) -> None:
    metadata_df = pd.read_csv(metadata_tsv, sep="\t", index_col=0)
    rna_df = pd.read_csv(rna_tsv, sep="\t", nrows=0)  # Just want the column names
    cell_type_cols = map_columns_to_cell_type(list(rna_df.columns)[1:], metadata_df)
    cell_type_files = create_cell_type_rna_files(rna_df, cell_type_cols, output_folder)

    with pd.read_csv(rna_tsv, sep="\t", skiprows=1, chunksize=CHUNK_SIZE, names=rna_df.columns) as reader:
        for chunk in reader:
            for cell_type, cols in cell_type_cols.items():
                chunk[cols].to_csv(cell_type_files[cell_type], index=False, header=False, mode='a')
        
if __name__ == "__main__":
    main()

