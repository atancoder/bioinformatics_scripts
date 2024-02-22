How to create cell type specific fragments from a scATAC folder (e.g http://catlas.org/catlas_downloads/humantissues/fragments)  


Create a new CellMetadata.tsv file, with a new column corresponding to fragment file
```
py attach_frag_file.py --metadata /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/Cell_metadata.tsv.gz --fragments /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/fragment_files --output cell_metadata_modified.tsv.gz
```

Takes a while. ~90 min
```
py create_fragments_per_cell_type.py --metadata cell_metadata_modified.tsv.gz --output_folder /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/cell_type_fragment_files > create_fragments_per_cell_type.log
```

Count fragments per cell
```
py count_num_fragments.py --cell_fragments /oak/stanford/groups/engreitz/Users/atan5133/data/CZI_scATAC/bingren_atlas/cell_type_fragment_files --output cell_frag_counts.tsv
```

Filter for cells with sufficient coverage
```
py filter_cell_types.py --cell_frag_counts cell_frag_counts.tsv --output cell_frag_counts_filtered.tsv
```