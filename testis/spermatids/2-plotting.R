# Plotting the trajectory of spermatocytes
# working directory is the same directory as the script

library(ComplexHeatmap)
library(tidyverse)
library(dynwrap)
library(dynplot)

source("0-common.R")

dataset <- read_rds("dataset.rds")
trajectory <- read_rds("trajectory.rds")
seu_oi <- read_rds("seu_oi.rds")
feature_importances <- read_rds("feature_importance.rds")

features_oi <- as.character(feature_importances$feature_id %>% head(5000))
# features_oi <- unique(c(features_oi, c("kmg")))

linearised <- linearise_cells(trajectory)
progressions <- linearised$progressions# %>% group_by(from, to) %>% sample_n(100, replace = TRUE) %>% arrange(cumpercentage)
expression_matrix <- t(as.matrix(seu_oi@assays$RNA@data[features_oi, progressions$cell_id]))
dim(expression_matrix)

# order features
clustering <- hclust(as.dist(dynutils::correlation_distance(t(expression_matrix))), method = "average")
clust <- as.dendrogram(clustering)

expression_matrix <- Matrix::t(dynutils::scale_quantile(expression_matrix))

feature_scores <- rank_features(expression_matrix, 0.9)

clust <- stats::reorder(clust, feature_scores, agglo.FUN = mean)

# Scatter + line ----------------------------------------------------------

colnames(seu_oi@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
plotdata_cells <- bind_cols(
  seu_oi@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell"),
  seu_oi@meta.data
)
colnames(plotdata_cells)

# project
layout <- trajectory %>% dynwrap::project_trajectory(seu_oi@reductions$umap@cell.embeddings, trajectory_projection_sd  = 2)
colnames(layout$dimred_segment_points) <- colnames(seu_oi@reductions$umap@cell.embeddings)

segment_info <- layout$dimred_segment_progressions %>%
  mutate(point_id = rownames(layout$dimred_segment_points)) %>% 
  mutate(edge_id = paste0(from, "->", to)) %>%
  arrange(edge_id, percentage) %>%
  left_join(layout$dimred_segment_points %>% as.data.frame() %>% rownames_to_column("point_id"), "point_id") %>%
  left_join(trajectory$milestone_network %>% select(-length), c("from", "to"))

plot <- ggplot(plotdata_cells, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = annotation)) +
  geom_path(aes(group = edge_id), data = segment_info, size = 2) +
  scale_cluster_color +
  theme_void() +
  theme(legend.position = "none")

ggsave(plot = plot, filename = "dimred.pdf")


# Smoothed expression -----------------------------------------------------

trajectory$milestone_network$waypoints <- pmap(
  trajectory$milestone_network,
  function(from, to, length, ...){
    tibble(
      percentage = seq(0, 1, length.out = length * 10), cell_id = paste0("W_", from, "_", to, "_", seq(length.out = length * 10)))
  }
)
waypoint_progressions <- trajectory$milestone_network %>% unnest(waypoints)
waypoint_milestone_percentages <- dynwrap::convert_progressions_to_milestone_percentages(waypoint_progressions$cell_id, trajectory$milestone_ids, trajectory$milestone_network, waypoint_progressions)
waypoint_milestone_percentages <- dplyr::rename(waypoint_milestone_percentages, waypoint_id = cell_id)

cell_ids <- trajectory$cell_ids[1:1000]
trajectory2 <- dynwrap::wrap_data(cell_ids = cell_ids) %>% add_trajectory(milestone_network = trajectory$milestone_network %>% select(from, to, directed, length), progressions = trajectory$progressions %>% filter(cell_id %in% cell_ids))
geodesic <- dynwrap::calculate_geodesic_distances(
  trajectory2,
  waypoint_milestone_percentages = waypoint_milestone_percentages
)

bw <- 1.

weights <- dnorm(geodesic, sd = bw)
weights <- weights / rowSums(weights)
weights <- weights[waypoint_progressions$cell_id, ]

# Annotation (celltypes)
celltype <- factor(seu_oi$annotation[colnames(weights)], levels = clusters)
celltype_matrix <- model.matrix(~0+celltype)
colnames(celltype_matrix) <- levels(celltype)

weights2 <- abs(geodesic) <= bw
weights2 <- weights2 / rowSums(weights2)
weights2 <- weights2[waypoint_progressions$cell_id, ]

celltype_smoothed <- weights2 %*% celltype_matrix
celltype_smoothed <- celltype_smoothed / rowSums(celltype_smoothed)

annotation_celltype <- anno_barplot(
  celltype_smoothed,
  gp = gpar(fill = cluster_colors[colnames(celltype_smoothed)], col = "#FFFFFF00"),
  bar_width = 1,
  height = unit(1, "cm"),
  axis = FALSE,
  axis_param = list(at = NULL)
)

# expression
expression_matrix <- t(as.matrix(seu_oi@assays$RNA@data[features_oi, trajectory2$cell_id]))

dim(expression_matrix)
dim(weights)

smoothed <- weights %*% expression_matrix

smoothed <- smoothed[, (apply(smoothed, 2, max) - apply(smoothed, 2, min)) > 0.5]

dim(smoothed)

scaled <- t(dynutils::scale_minmax(smoothed))
feature_scores <- rank_features(scaled, 0.9)
scaled <- scaled[names(feature_scores)[feature_scores %>% order()],]

# cluster/order
clustering <- hclust(as.dist(dynutils::correlation_distance(scaled)), method = "average")
clust <- as.dendrogram(clustering)
clust <- stats::reorder(clust, feature_scores[rownames(scaled)], agglo.FUN = mean)

# Library size
library <- Matrix::rowSums(dataset$counts)[trajectory2$cell_id]
library_smoothed <- weights %*% t(t(library))

max_lib <- ceiling(max(library_smoothed) / 1000) * 1000
annotation_library <- anno_lines(library_smoothed, axis_param = list(at = c(0, max_lib)), ylim = c(0, max_lib))

dim(t(library))
dim(weights)

## Annotation
heatmap_legend_param <- list(
  direction = "horizontal",
  border = "#333333",
  at = c(0, 1),
  labels = c("low", "high"),
  legend_width = unit(3, "cm"),
  title_position = "topcenter"
)

symbols_oi <- feature_importances %>% 
  filter(!(str_detect(feature_id, "asRNA:") | str_detect(feature_id, "lncRNA:") | str_detect(feature_id, "CG"))) %>% 
  left_join(tibble(
    feature_id = names(feature_scores),
    bin = feature_scores %>% cut(breaks = 5),
    # bin = cutree(clustering, 10)
  )) %>% 
  group_by(bin) %>% 
  top_n(2, wt = importance) %>% 
  pull(feature_id)
symbols_oi <- c(
  symbols_oi,
  "kmg", "Rbp4", "fzo", "can", "sa", "mst87F", "kl-3", "kl-5", "eya",
  NULL
)
left_annotation <- rowAnnotation(
  gene = anno_mark(at = match(symbols_oi, rownames(scaled)), labels = symbols_oi, side  = "left")
)

top_annotation <- columnAnnotation(
  `# UMIs` = annotation_library,
  annotation = annotation_celltype
)


heatmap <- Heatmap(
  # scaled,
  scaled,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  top_annotation = top_annotation,
  # top_annotation = column_annotation,
  left_annotation = left_annotation,
  col = rev(RColorBrewer::brewer.pal(10, "RdBu")[2:9]),
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = heatmap_legend_param,
  show_row_dend = FALSE,
  name = "Expression",
  use_raster = TRUE,
  raster_quality = 5
)
heatmap_list <- ComplexHeatmap::HeatmapList()
heatmap_list <- ComplexHeatmap::add_heatmap(heatmap, heatmap_list)

annotation_legends = c(
  # column_annotation_legends
)

heatmap <- ComplexHeatmap::make_layout(
  heatmap_list,
  annotation_legend_list = annotation_legends,
  merge_legend = TRUE,
  heatmap_legend_side = "top",
  annotation_legend_side = "top"
)
# heatmap

pdf("heatmap.pdf", width = 6, height = 6);draw(heatmap);dev.off()

dim(smoothed)

