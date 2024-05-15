import concurrent.futures
import glob
import gzip
import os
import re
import subprocess
import time

import click
import pandas as pd

NUM_THREADS = 32


def find_replicate_number(filename):
    pattern = r"_rep(\d+)_"
    # Search the path for the pattern
    match = re.search(pattern, filename)

    # Check if a match was found and print the result
    if match:
        captured = int(match.group(1))
        return captured
    else:
        raise Exception(f"Cannot infer replicate number in {filename}")


def split_sample(frag_file, metadata_df, chromosomes, output_folder):
    sample_name = metadata_df[metadata_df["frag_file"] == frag_file].iloc[0]["sample"]
    cell_type_files = {}
    fragments_written = 0
    # replicate = find_replicate_number(frag_file)
    with gzip.open(frag_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, start, end, barcode, read_count = line.split("\t")
            if chrom not in chromosomes:
                continue
            cell_id = f"{sample_name}#{barcode}"
            if cell_id in metadata_df.index:
                row = metadata_df.loc[cell_id]
                cell_type = row["cell_type"]
                if cell_type not in cell_type_files:
                    cell_type_files[cell_type] = gzip.open(
                        os.path.join(
                            output_folder, f"{sample_name}_{cell_type}.fragments.bed.gz"
                        ),
                        "wt",
                    )
                barcode = f"{sample_name}.{barcode}".replace('-', '.')
                line = "\t".join([chrom, start, end, barcode, read_count])
                cell_type_files[cell_type].write(line)
                fragments_written += 1

    for _, f in cell_type_files.items():
        f.close()

    return fragments_written


def split_samples(frag_files, metadata_df, chromosomes, output_folder):
    """
    split a sample fragment file into multiple cell type fragment files
    """
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = {
            executor.submit(split_sample, frag_file, metadata_df, chromosomes, output_folder)
            for frag_file in frag_files
        }

        for _, future in enumerate(concurrent.futures.as_completed(futures)):
            future.result()


def merge_frag_files(frag_files, cell_type_filename):
    frag_files_str = " ".join(frag_files)
    cmd = f"cat {frag_files_str} > {cell_type_filename} && rm {frag_files_str}"
    subprocess.run(cmd, shell=True, check=True)


def merge_cell_types(output_folder, metadata_df):
    """
    Merge fragment files with the same cell type
    """
    frag_files_to_merge = []
    for cell_type in metadata_df["cell_type"].drop_duplicates():
        frag_files = glob.glob(
            os.path.join(output_folder, f"*_{cell_type}.fragments.bed.gz")
        )
        if not frag_files:
            continue
        cell_type_fragment = os.path.join(
            output_folder, f"{cell_type}.fragments.bed.gz"
        )
        frag_files_to_merge.append([frag_files, cell_type_fragment])

    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = {
            executor.submit(merge_frag_files, frag_files, cell_type_filename)
            for frag_files, cell_type_filename in frag_files_to_merge
        }

        for _, future in enumerate(concurrent.futures.as_completed(futures)):
            future.result()


@click.command
@click.option("--metadata", "metadata_tsv", required=True)
@click.option("--chromosomes", required=True)
@click.option("--output_folder", required=True)
def main(metadata_tsv, chromosomes, output_folder) -> None:
    metadata_df = pd.read_csv(metadata_tsv, sep="\t", index_col=0)
    chromosomes = set(pd.read_csv(chromosomes, sep='\t', names=['chr', 'size'])['chr'])
    frag_files = list(
        metadata_df[~metadata_df["frag_file"].isna()]["frag_file"].drop_duplicates()
    )

    start_time = time.time()
    split_samples(frag_files, metadata_df, chromosomes, output_folder)
    print(f"Split samples in {time.time() - start_time} s")

    start_time = time.time()
    merge_cell_types(output_folder, metadata_df)
    print(f"Merge cell types in {time.time() - start_time} s")


if __name__ == "__main__":
    main()
