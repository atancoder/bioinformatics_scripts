import click
import pandas as pd
import glob
import time
import os
import re
import gzip
import concurrent.futures
import subprocess

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


def split_tissue(frag_file, metadata_df, output_folder):
    tissue_name = metadata_df[metadata_df["frag_file"] == frag_file].iloc[0]["tissue"]
    cell_type_files = {}
    fragments_written = 0
    replicate = find_replicate_number(frag_file)
    with gzip.open(frag_file, "rt") as f:
        for line in f:
            barcode = line.split("\t")[3]
            cell_id = f"{tissue_name}_{replicate}+{barcode}"
            if cell_id in metadata_df.index:
                row = metadata_df.loc[cell_id]
                cell_type = row["cell_type"]
                if cell_type not in cell_type_files:
                    cell_type_files[cell_type] = gzip.open(
                        os.path.join(
                            output_folder, f"{tissue_name}_{cell_type}.fragments.bed.gz"
                        ),
                        "wt",
                    )
                cell_type_files[cell_type].write(line)
                fragments_written += 1

    for _, f in cell_type_files.items():
        f.close()

    return fragments_written


def split_tissues(frag_files, metadata_df, output_folder):
    """
    split a tissue fragment file into multiple cell type fragment files
    """
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = {
            executor.submit(split_tissue, frag_file, metadata_df, output_folder)
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
@click.option("--output_folder", required=True)
def main(metadata_tsv, output_folder) -> None:
    metadata_df = pd.read_csv(metadata_tsv, sep="\t", index_col=0)
    frag_files = list(
        metadata_df[~metadata_df["frag_file"].isna()]["frag_file"].drop_duplicates()
    )

    start_time = time.time()
    split_tissues(frag_files, metadata_df, output_folder)
    print(f"Split tissues in {time.time() - start_time}")

    start_time = time.time()
    merge_cell_types(output_folder, metadata_df)
    print(f"Merge cell types in {time.time() - start_time}")


if __name__ == "__main__":
    main()
