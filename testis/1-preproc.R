library(loomR)
library(Seurat)
library(hdf5r)
library(tidyverse)

setwd("testis")

# download the v1 loom file from SCOPE
lfile <- connect(("r_fca_biohub_testis_10x.loom"), mode = "r", skip.validate = T)

counts <- lfile[["matrix"]][,]

# metadata of cells
col_attrs <- list.datasets(lfile[["col_attrs"]]) %>% map(function(attr) {
  lfile[["col_attrs"]][[attr]][]
})
names(col_attrs) <- list.datasets(lfile[["col_attrs"]])
metadata <- col_attrs %>% keep(is.vector) %>% as_tibble()
metadata$cell_id <- metadata$CellID

rownames(counts) <- metadata$cell_id

# metadata of genes
row_attrs <- list.datasets(lfile[["row_attrs"]]) %>% map(function(attr) {
  lfile[["row_attrs"]][[attr]][]
})
names(row_attrs) <- list.datasets(lfile[["row_attrs"]])
feature_metadata <- row_attrs %>% keep(is.vector) %>% as_tibble()
feature_metadata$feature_id <- feature_metadata$Gene

colnames(counts) <- feature_metadata$feature_id

# Create seurat
seu <- Seurat::CreateSeuratObject(Matrix::t(counts), meta.data = column_to_rownames(as.data.frame(metadata), "cell_id"))

seu <- seu %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

embedding <- lfile[["col_attrs"]][["Embedding"]][]
colnames(embedding) <- c("UMAP_1", "UMAP_2")
rownames(embedding) <- colnames(seu)
umap <- CreateDimReducObject(as.matrix(embedding), key = "umap_", assay = "RNA")
seu@reductions$umap <- umap

# add desired clustering and annotation

write_rds(seu, "testis.rds")
