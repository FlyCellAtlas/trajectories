library(Seurat)
source("R/experiment.R")
library(tidyverse)

devtools::load_all("~/thesis/projects/dynverse/dynplot2/")
devtools::load_all("~/thesis/projects/dynverse/dynplot/")
devtools::load_all("~/thesis/projects/dynverse/libraries/dynfeature")

library(dynwrap)

experiment_testis <- Experiment$new("testis")
experiment <- experiment_testis$subexperiment("cysts")
experiment <- experiment_testis$subexperiment("cysts_complex")

seu <- read_rds(experiment_testis$output_file("testis.rds"))
DimPlot(seu, group.by = "annotation", label = TRUE)

clusters <- c(
  "cyst stem cell",
  "early cyst cells",
  "Spermatocyte Cyst Cell B",
  "Spermatocyte Cyst Cell A",
  # "cyst cell intermediate",
  "cyst cell branch",
  # "late cyst cell branch a",
  # "late cyst cell branch b",
  "terminal epithelial cell of testis",
  NULL
)

seu_oi <- seu[,seu$annotation %in% clusters]

seu_oi$annotation_upper <- c(
  "cyst stem cell" = "cyst stem cell",
  "early cyst cells" = "early cyst cells",
  "Spermatocyte Cyst Cell B" = "spermatocyte cyst cell",
  "Spermatocyte Cyst Cell A" = "spermatocyte cyst cell",
  "cyst cell branch" =  "spermatocyte cyst cell",
  "cyst cell intermediate" = "cyst cell intermediate",
  "late cyst cell branch a" = "late cyst cell",
  "late cyst cell branch b" = "late cyst cell",
  "terminal epithelial cell of testis" = "??"
)[seu_oi$annotation]

plot_dimred <- DimPlot(seu_oi, reduction = "umap", group.by = "annotation") + scale_cluster_color
plot_dimred
save_plot(plot_dimred, experiment$results_file("dimred_original"), width = 8, height = 5)

# dimplot <- DimPlot(seu_oi, group.by = "annotation") + scale_cluster_color

# seu_oi <- seu_oi[, seu_oi$seurat_clusters %in% c(0, 6)]

seu_oi <- seu_oi %>% NormalizeData() %>% FindVariableFeatures()
seu_oi <- seu_oi %>% ScaleData()

seu_oi <- seu_oi %>% RunPCA()
seu_oi@reductions$pca@feature.loadings %>% dim

DimPlot(seu_oi, reduction = "pca", group.by = "annotation", label = TRUE)


seu_oi <- seu_oi %>% RunUMAP(dims = 1:50)
seu_oi <- seu_oi %>% FindNeighbors() %>% FindClusters(resolution = 0.4)
seu_oi <- seu_oi[,seu_oi$seurat_clusters != 7]

DimPlot(seu_oi, reduction = "pca", group.by = "seurat_clusters", label = TRUE)

seu_oi$log_nCount_RNA <- log10(seu_oi$nCount_RNA)
FeaturePlot(seu_oi, reduction = "pca", features = "log_nCount_RNA", label = TRUE)
DimPlot(seu_oi, reduction = "pca", group.by = "annotation", label = TRUE) + scale_cluster_color

plot_dimred <- DimPlot(seu_oi, reduction = "umap", group.by = "annotation", label = FALSE) + scale_cluster_color
plot_dimred
save_plot(plot_dimred, experiment$results_file("dimred"), width = 8, height = 5)

FeaturePlot(seu_oi, reduction = "umap", features = "nCount_RNA")

umap <- seu_oi@reductions$umap@cell.embeddings
pca <- seu_oi@reductions$pca@cell.embeddings

dataset <- wrap_expression(
  expression = t(as.matrix(seu_oi@assays$RNA@data[seu_oi@assays$RNA@var.features,])),
  counts = t(as.matrix(seu_oi@assays$RNA@counts[seu_oi@assays$RNA@var.features,]))
)
mds <- dyndimred::dimred_landmark_mds(dataset$expression, ndim = 3)
dataset <- dataset %>% add_prior_information(
  # dimred = umap,
  # dimred = mds[,],
  dimred = as.matrix(seu_oi@reductions$pca@cell.embeddings[,1:20]),
  groups_id = seu_oi$seurat_clusters %>% enframe("cell_id", "group_id")
  # groups_id = seu_oi$annotation %>% enframe("cell_id", "group_id")
  # groups_id = seu_oi$annotation_upper %>% enframe("cell_id", "group_id")
)
dataset$cell_info$groups_id <- dataset$prior_information$groups_id$group_id

dataset <- dataset %>% add_dimred(umap)
dataset <- dataset %>% add_dimred(pca)
dataset <- dataset %>% add_dimred(mds)
dynplot_dimred(dataset) + 
  geom_cell_point(aes(color = groups_id)) + 
  geom_cell_contour(aes(group = groups_id, fill = groups_id)) +
  geom_cell_contour_label(aes(group = groups_id, label = groups_id, color = groups_id)) +
  # scale_cluster_fill +
  # scale_cluster_color +
  theme(legend.position = "none")

# method <- SCORPIUS::ti_scorpius()
method <- tislingshot::ti_slingshot()

##
library(slingshot)
rd <- dataset$prior_information$dimred
labels <- dataset$prior_information$groups_id %>% deframe()
sds <- slingshot::slingshot(
  rd,
  labels,
  start.clus = "cyst stem cell"
)
start_cell <- apply(slingshot::slingPseudotime(sds), 1, min) %>% sort() %>% head(1) %>% names()
start.clus <- labels[[start_cell]]

lineages <- slingLineages(sds)
lineage_ctrl <- slingParams(sds)
cluster_network <- lineages %>%
  map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
  unique() %>%
  mutate(
    length = lineage_ctrl$dist[cbind(from, to)],
    directed = TRUE
  )

# collect dimred
dimred <- reducedDim(sds)

# collect clusters
cluster <- slingClusterLabels(sds)

# collect progressions
adj <- slingAdjacency(sds)
lin_assign <- apply(slingCurveWeights(sds), 1, which.max)

progressions <- map_df(seq_along(lineages), function(l) {
  ind <- lin_assign == l
  lin <- lineages[[l]]
  pst.full <- slingPseudotime(sds, na = FALSE)[,l]
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

#   ____________________________________________________________________________
#   Save output                                                             ####
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

##

dynplot::plot_topology(trajectory)# + scale_cluster_fill

colnames(umap) <- paste0("comp_", seq(ncol(umap)))
colnames(pca) <- paste0("comp_", seq(ncol(pca)))


# use umap
traj_dimred_umap <- trajectory %>% dynwrap::project_trajectory(umap)
trajectory[names(traj_dimred_umap)] <- traj_dimred_umap
trajectory$dimred <- umap

# use pca
# traj_dimred_pca <- trajectory %>% dynwrap::project_trajectory(pca)
# trajectory[names(traj_dimred_pca)] <- traj_dimred_pca
# trajectory$dimred <- pca

devtools::load_all("~/thesis/projects/dynverse/dynplot2/")
devtools::load_all("~/thesis/projects/dynverse/dynplot/")

plot_dimred(trajectory, hex_cells = FALSE)
plot_dimred(trajectory, hex_cells = FALSE, color_cells = "pseudotime")
dynplot::plot_topology(trajectory)

trajectory$milestone_network

feature_importances <- dynfeature::calculate_overall_feature_importance(trajectory, dataset)
# milestone_feature_importances <- dynfeature::calculate_milestone_feature_importance(trajectory, dataset)


write_rds(dataset, experiment$output_file("dataset.rds"))
write_rds(trajectory, experiment$output_file("trajectory.rds"))
write_rds(seu_oi, experiment$output_file("seu_oi.rds"))
write_rds(umap, experiment$output_file("umap.rds"))
write_rds(pca, experiment$output_file("pca.rds"))
write_rds(feature_importances, experiment$output_file("feature_importance.rds"))

