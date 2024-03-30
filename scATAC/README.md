How to create cell type specific fragments from a scATAC folder (e.g http://catlas.org/catlas_downloads/humantissues/fragments)  


Create a new CellMetadata.tsv file, with a new column corresponding to fragment file and other processing done to names
```
py preprocess_metadata.py --metadata /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/Cell_metadata.tsv.gz --fragments /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/fragment_files --output cell_metadata_modified.tsv.gz
```

Takes a while. ~80 min
```
py create_fragments_per_cell_type.py --metadata cell_metadata_modified.tsv.gz --output_folder /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/cell_type_fragment_files --chromososomes /oak/stanford/groups/engreitz/Users/atan5133/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv
```

Count fragments per cell (for relevant chromosomes)
```
py count_num_fragments.py --cell_fragments /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/cell_type_fragment_files --output cell_frag_counts.tsv  --chrom_bed /oak/stanford/groups/engreitz/Users/atan5133/anthonys_scripts/scATAC/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
```

Filter for cells with sufficient coverage
```
py filter_cell_types.py --cell_frag_counts cell_frag_counts.tsv --output cell_frag_counts_filtered.tsv
```

Convert files to tagAlign
Modify the following config to adjust default file locations if needed: `snakefile_config.yml`
```
snakemake --snakefile Snakefile_tagAlign --profile slurm
```

Create biosamples.config for ABC
```
py gen_biosamples.py --tagAlign_folder /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/cell_type_tagAlign_files --output biosamples_config.tsv
```
