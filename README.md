This contains the scripts to perform the trajectory inference on _spermatocytes_, _spermatids_ and _early cyst cells_.

### Dependencies

```
# CRAN:

loomR
Seurat
hdf5r
tidyverse
Seurat
dynwrap
Seurat

# Bioconductor:

slingshot
ComplexHeatmap
```

### Running the code

- Download the "relaxed" loom file from SCOPE (https://scope.aertslab.org/#/FlyCellAtlas) named `r_fca_biohub_testis_10x.loom` and put it in the testis directory
- Run the scripts in the testis directory in order
- Run the scripts in spermatocytes, spermatids or cysts directories