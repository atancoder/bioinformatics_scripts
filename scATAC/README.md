How to create cell type specific fragments from a scATAC folder (e.g http://catlas.org/catlas_downloads/humantissues/fragments)  


Create a new CellMetadata.tsv file, with a new column corresponding to fragment file and other processing done to names
```
py preprocess_metadata.py --metadata /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/Cell_metadata.tsv.gz --fragments /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/fragment_files --output catlas_cell_metadata.tsv.gz
```

Takes a while. ~2 hours with 16 cpus
```
py create_fragments_per_cell_type.py --metadata catlas_cell_metadata.tsv.gz --output_folder /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/cell_type_fragment_files --chromosomes GRCh38_EBV.no_alt.chrom.sizes.tsv
```

Count fragments per cell (for relevant chromosomes)
```
py count_num_fragments.py --cell_fragments /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/cell_type_fragment_files --output catlas_cell_frag_counts.tsv
```

Filter for cells with sufficient coverage
```
py filter_cell_types.py --cell_frag_counts catlas_cell_frag_counts.tsv --output catlas_cell_frag_counts_filtered.tsv
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

## For Stitziel Multiome
```
py preprocess_metadata.py --metadata /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/meta/coronary_multiome_meta.tsv --fragments /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/fragments --output /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/stitziel_multiome_meta.tsv.gz
```

```
py create_fragments_per_cell_type.py --metadata /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/stitziel_multiome_meta.tsv.gz --output_folder /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/cell_type_fragment_files --chromosomes GRCh38_EBV.no_alt.chrom.sizes.tsv
```

To split an aggregated RNA matrix into cell type specific RNA matrices
```
py split_rna_matrix.py --rna /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/counts/coronary_multiome_RNA.tsv --metadata /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/stitziel_multiome_meta.tsv.gz --output_folder /oak/stanford/groups/engreitz/Users/atan5133/data/stitziel_multiome/cell_type_rna_matrices
```


