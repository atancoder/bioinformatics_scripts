import click
import os
import concurrent.futures
import subprocess
import pandas as pd

NUM_THREADS = 32

def count_fragment_file(filename, chrom_bed):
    cmd = f"bedtools intersect -nonamecheck -a {filename} -b {chrom_bed} | wc -l"
    return int(subprocess.check_output(cmd, shell=True, text=True).split(' ')[0])

@click.command
@click.option("--cell_fragments", "cell_fragments_folder", required=True)
@click.option("--chrom_bed", required=True)
@click.option("--output", required=True)
def main(cell_fragments_folder, chrom_bed, output):
    files = os.listdir(cell_fragments_folder)
    frag_files = [os.path.join(cell_fragments_folder, f) for f in files if "fragments" in f]
    print(f"Found {len(frag_files)} fragment files")
    cell_type_fragments = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = {
            executor.submit(count_fragment_file, chrom_bed, frag_file)
            for frag_file in frag_files
        }

        for idx, future in enumerate(concurrent.futures.as_completed(futures)):
            num_fragments = future.result() 
            cell_type_frag_file = os.path.basename(frag_files[idx])
            cell_type = cell_type_frag_file.split('.')[0]
            cell_type_fragments.append({"cell_type": cell_type, "file": frag_files[idx], "num_fragments": num_fragments})
    
    pd.DataFrame(cell_type_fragments).to_csv(output, sep='\t', index=False)

if __name__ == "__main__":
    main()