# Inferring the trajectory for spermatocytes
# working directory is the same directory as the script

library(Seurat)
library(tidyverse)
library(slingshot)
library(dynwrap)

setwd(rprojroot::find_root(rprojroot::is_git_root))
setwd('testis/spermatids/')

source("0-common.R")

seu <- read_rds("../testis.rds")
DimPlot(seu, group.by = "annotation", label = TRUE)

clusters <- c(
  "early elongation stage spermatid",
  "early-mid elongation-stage spermatid",
  "mid-late elongation-stage spermatid"
)

seu_oi <- seu[,seu$annotation %in% clusters]

seu_oi@reductions$umap_orig <- seu_oi@reductions$umap

plot_dimred <- DimPlot(seu_oi, group.by = "annotation") + scale_cluster_color
plot_dimred

seu_oi <- seu_oi %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:50)
seu_oi <- seu_oi %>% FindNeighbors() %>% FindClusters(resolution = 0.4)

DimPlot(seu_oi, reduction = "pca", group.by = "seurat_clusters", label = TRUE)
DimPlot(seu_oi, reduction = "pca", group.by = "annotation", label = TRUE)
plot_dimred <- DimPlot(seu_oi, reduction = "umap", group.by = "annotation", label = TRUE) + scale_cluster_color
plot_dimred

DimPlot(seu_oi, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(seu_oi, reduction = "umap", group.by = "annotation", label = TRUE)
FeaturePlot(seu_oi, reduction = "umap", features = "nCount_RNA")

plot_dimred <- DimPlot(seu_oi, reduction = "umap", group.by = "annotation", label = TRUE) + scale_cluster_color
plot_dimred

plot_features <- FeaturePlot(seu_oi, features = c("soti", "dac", "zormin"), ncol = 3)
plot_features

umap <- seu_oi@reductions$umap@cell.embeddings
pca <- seu_oi@reductions$pca@cell.embeddings

# create dynwrap dataset
dataset <- wrap_expression(
  expression = t(as.matrix(seu_oi@assays$RNA@data[seu_oi@assays$RNA@var.features,])),
  counts = t(as.matrix(seu_oi@assays$RNA@counts[seu_oi@assays$RNA@var.features,]))
)
mds <- dyndimred::dimred_landmark_mds(dataset$expression, ndim = 2)
dataset <- dataset %>% add_prior_information(
  dimred = as.matrix(seu_oi@reductions$pca@cell.embeddings[,1:10]), # use the first 10 PCA dimensions
  groups_id = seu_oi$annotation %>% enframe("cell_id", "group_id") # use the annotation for clustering
)
dataset$cell_info$groups_id <- dataset$prior_information$groups_id$group_id

dataset <- dataset %>% add_dimred(mds)

# Run slingshot -----------------------------------------------------------
rd <- dataset$prior_information$dimred
labels <- dataset$prior_information$groups_id %>% deframe()
sds <- slingshot::slingshot(
  rd,
  labels,
  start.clus = "early elongation stage spermatid"
)
start_cell <- apply(slingshot::slingPseudotime(sds), 1, min) %>% sort() %>% head(1) %>% names()
start.clus <- labels[[start_cell]]

lineages <- slingshot::slingLineages(sds)
cluster_network <- lineages %>%
  map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
  unique() %>%
  mutate(
    length = map2_dbl(from, to, ~as.numeric(igraph::as_adjacency_matrix(sds@metadata$mst, attr = "weight")[.x, .y])),
    directed = TRUE
  )

# collect clusters
cluster <- slingshot::slingClusterLabels(sds)

# collect progressions
lin_assign <- apply(slingshot::slingCurveWeights(sds), 1, which.max)

progressions <- map_df(seq_along(lineages), function(l) {
  ind <- lin_assign == l
  lin <- lineages[[l]]
  pst.full <- slingshot::slingPseudotime(sds, na = FALSE)[,l]
  pst <- pst.full[ind]
  means <- sapply(lin, function(clID){
    stats::weighted.mean(pst.full, cluster[,clID])
  })
  non_ends <- means[-c(1,length(means))]
  edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
  from.l <- lineages[[l]][edgeID.l]
  to.l <- lineages[[l]][edgeID.l + 1]
  m.from <- means[from.l]
  m.to <- means[to.l]
  
  pct <- (pst - m.from) / (m.to - m.from)
  pct[pct < 0] <- 0
  pct[pct > 1] <- 1
  
  tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
})



# Create dynwrap object ---------------------------------------------------

trajectory <-
  dynwrap::wrap_data(
    cell_ids = rownames(rd)
  ) %>%
  dynwrap::add_trajectory(
    milestone_network = cluster_network,
    progressions = progressions
  ) %>%
  dynwrap::add_dimred(
    dimred = dataset$prior_information$dimred
  )

dynplot::plot_topology(trajectory)

# simplify
trajectory <- trajectory %>% simplify_trajectory()

# determine differential genes
feature_importances <- dynfeature::calculate_overall_feature_importance(trajectory, dataset)

# Store -------------------------------------------------------------------
write_rds(dataset, "dataset.rds")
write_rds(trajectory, "trajectory.rds")
write_rds(seu_oi, "seu_oi.rds")
write_rds(feature_importances, "feature_importance.rds")
