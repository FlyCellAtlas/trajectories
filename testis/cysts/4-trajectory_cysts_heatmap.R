library(ComplexHeatmap)
library(Seurat)
library(tidyverse)
source("R/experiment.R")
source("code/spermatogenesis_v2/4-common.R")


devtools::load_all("~/thesis/projects/dynverse/dynplot2/")
devtools::load_all("~/thesis/projects/dynverse/dynplot/")
devtools::load_all("~/thesis/projects/dynverse/libraries/dynfeature")

library(dynwrap)

experiment_testis <- Experiment$new("testis")
experiment <- experiment_testis$subexperiment("cysts")

dataset <- read_rds(experiment$output_file("dataset.rds"))
trajectory <- read_rds(experiment$output_file("trajectory.rds"))
seu_oi <- read_rds(experiment$output_file("seu_oi.rds"))
annot_mapping <- read_rds(experiment_testis$output_file("annot_mapping.rds"))
seu_oi$annotation <- annot_maping[colnames(seu_oi)]
umap <- read_rds(experiment$output_file("umap.rds"))
pca <- read_rds(experiment$output_file("pca.rds"))
feature_importances <- read_rds(experiment$output_file("feature_importance.rds"))

trajectory <- trajectory %>% simplify_trajectory()
trajectory <- trajectory %>% add_root(root_cell_id = colnames(seu_oi)[seu_oi$annotation == "cyst stem cell"][1])

plot_topology(trajectory)

trajectory$milestone_network


# Smoothed expression -----------------------------------------------------

# Waypoints
milestone_network <- trajectory$milestone_network
milestone_network$waypoints <- pmap(
  milestone_network,
  function(from, to, length, ...){
    tibble(
      percentage = seq(0, 1, length.out = length * 10), cell_id = paste0("W_", from, "_", to, "_", seq(length.out = length * 10)))
  }
)
waypoint_progressions <- milestone_network %>% unnest(waypoints)
waypoint_progressions$edge_ix <- group_indices(waypoint_progressions %>% group_by(from, to))
waypoint_milestone_percentages <- dynwrap::convert_progressions_to_milestone_percentages(waypoint_progressions$cell_id, milestone_ids, trajectory$milestone_network, waypoint_progressions)
waypoint_milestone_percentages <- rename(waypoint_milestone_percentages, waypoint_id = cell_id)

column_split <- waypoint_progressions$edge_ix

cell_ids <- trajectory$cell_ids#[1:1000]
trajectory2 <- dynwrap::wrap_data(cell_ids = cell_ids) %>% add_trajectory(milestone_network = milestone_network %>% select(from, to, directed, length), progressions = trajectory$progressions %>% filter(cell_id %in% cell_ids))
geodesic <- dynwrap::calculate_geodesic_distances(
  trajectory2,
  waypoint_milestone_percentages = waypoint_milestone_percentages
)

bw <- mean(geodesic) / 10

weights <- dnorm(geodesic, sd = bw)
weights <- weights / rowSums(weights)
weights <- weights[waypoint_progressions$cell_id, ]

# Annotation (celltypes)
celltype <- factor(seu_oi$annotation[colnames(weights)])
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

# Expression
features_oi <- as.character(feature_importances$feature_id %>% head(1000)) %>% c("stg") %>% unique()
"stg" %in% features_oi
expression_matrix <- t(as.matrix(seu_oi@assays$RNA@data[features_oi, trajectory2$cell_id]))

dim(expression_matrix)
dim(weights)

smoothed <- weights %*% expression_matrix

smoothed <- smoothed[, (apply(smoothed, 2, max) - apply(smoothed, 2, min)) > 0.5]

scaled <- t(dynutils::scale_minmax(smoothed))
feature_scores <- rank_features(scaled, 0.99)
scaled <- scaled[names(feature_scores)[feature_scores %>% order()],]



# pheatmap::pheatmap(t(scaled)[,order(feature_scores)], cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = FALSE)

# Gene clustering and ordering
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

# gene marks
max_edge <- aggregate(smoothed, list(waypoint_progressions$edge_ix), max)
ratio <- (apply(max_edge, 2, median) / apply(max_edge, 2, max)) < 0.5
genes_specific <- names(ratio)[ratio]
symbols_oi <- feature_importances %>% 
  filter(
    !(str_detect(feature_id, "asRNA:") | str_detect(feature_id, "lncRNA:") | str_detect(feature_id, "CG")) &
      feature_id %in% genes_specific
  ) %>% 
  left_join(tibble(
    feature_id = names(feature_scores),
    bin = feature_scores %>% cut(breaks = 10),
    # bin = cutree(clustering, 10)
  )) %>% 
  group_by(bin) %>% 
  top_n(2, wt = importance) %>% 
  pull(feature_id)
symbols_oi <- c(
  symbols_oi,
  "stg",
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
  column_split = factor(column_split, levels = unique(column_split)),
  left_annotation = left_annotation,
  col = rev(RColorBrewer::brewer.pal(10, "RdBu")[2:9]),
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = heatmap_legend_param,
  show_row_dend = FALSE,
  name = "Expression",
  column_title = NULL,
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

pdf(experiment$results_file("heatmap.pdf"), width = 6, height = 6);draw(heatmap);dev.off()

dim(scaled)


# Specific genes ----------------------------------------------------------


trajectory$milestone_network
waypoints_oi <- waypoint_progressions %>%
  filter(edge_ix %in% c(1, 2, 3)) %>%
  group_by(edge_ix) %>%
  arrange(percentage) %>%
  slice(-1) %>%
  pull(cell_id)
waypoint_to_edge_ix <- waypoint_progressions %>% select(cell_id, edge_ix) %>% deframe()
branch_expression <- aggregate((smoothed)[waypoints_oi,], list(edge = waypoint_to_edge_ix[waypoints_oi]), max) %>% pivot_longer(
  cols = colnames(smoothed),
  values_to = "expression",
  names_to = c("gene"),
)

start_edge <- 3
branch_specific <- branch_expression %>% 
  group_by(gene) %>% 
  summarize(
    stem_expression = expression[edge == start_edge],
    specificity = max(expression) - min(expression[edge != start_edge]),
    upregulation = max(expression) - stem_expression,
    edge_ix = edge[which.max(expression)],
    NULL
    ) %>% 
  ungroup()

edge_labelling <- list(
  `1`="Cyst cell branch",
  `2`="Spermatocyte Cyst Cell A"
)

branch_scores <- branch_specific %>% 
  filter(edge_ix != start_edge) %>% 
  arrange(desc(specificity * upregulation)) %>% 
  mutate(branch = as.character(edge_labelling[edge_ix])) %>% 
  select(gene, stem_expression, upregulation, specificity, branch)

writexl::write_xlsx(list(scores = branch_scores), "cyst_branch_genes.xlsx")



plot(smoothed[, "Trim9"])





# Scatter + line ----------------------------------------------------------


plotdata_cells <- bind_cols(
  seu_oi@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell"),
  seu_oi@meta.data
)

layout <- trajectory %>% dynwrap::project_trajectory(seu_oi@reductions$umap@cell.embeddings, trajectory_projection_sd  = 2)

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

save_plot(plot, experiment$results_file("scatter"), width = 6, height = 6)

plot +
  theme(legend.position = "right")
