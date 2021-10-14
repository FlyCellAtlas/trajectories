library(Seurat)
source("R/experiment.R")
library(tidyverse)

devtools::load_all("~/thesis/projects/dynverse/dynplot2/")
devtools::load_all("~/thesis/projects/dynverse/dynplot/")
devtools::load_all("~/thesis/projects/dynverse/libraries/dynfeature")

cluster_colors <- c(
  "cyst stem cell" = "#d53e4f",
  "early cyst cells" = "#fc8d59",
  "Spermatocyte Cyst Cell A" = "#fee08b",
  "Spermatocyte Cyst Cell B" = "#fdc116",
  "cyst cell branch" = "#e6f598",
  "cyst cell intermediate" = "#99d594",
  "late cyst cell branch a" = "#3288bd",
  "late cyst cell branch b" = "#81bbdd"
)


library(dynwrap)

experiment_testis <- Experiment$new("testis")
experiment <- experiment_testis$subexperiment("cysts")



dataset <- read_rds(experiment$output_file("dataset.rds"))
trajectory <- read_rds(experiment$output_file("trajectory.rds"))
seu_oi <- read_rds(experiment$output_file("seu_oi.rds"))
umap <- read_rds(experiment$output_file("umap.rds"))
pca <- read_rds(experiment$output_file("pca.rds"))
feature_importances <- read_rds(experiment$output_file("feature_importance.rds"))





library(ComplexHeatmap)
# devtools::install_github("jokergoo/ComplexHeatmap@c7865b7e791c3dd4d7bdd6ebb53f59567bd04111")

features_oi <- as.character(feature_importances$feature_id %>% head(100))
# features_oi <- milestone_feature_importances %>% 
#   group_by(milestone_id) %>% 
#   slice(1:10) %>% 
#   pull(feature_id) %>% 
#   as.character() %>% 
#   unique()

linearised <- linearise_cells(trajectory)
progressions <- linearised$progressions# %>% group_by(from, to) %>% sample_n(100, replace = TRUE) %>% arrange(cumpercentage)
expression_matrix <- as.matrix(dataset$expression[progressions$cell_id, features_oi])

# order features
clust <- hclust(as.dist(dynutils::correlation_distance(t(expression_matrix))), method = "average")
clust <- as.dendrogram(clust)

expression_matrix <- Matrix::t(dynutils::scale_quantile(expression_matrix))

feature_scores <- dynplot:::rank_features_quantile(expression_matrix, 0.8)

clust <- stats::reorder(clust, feature_scores, agglo.FUN = mean)

# column annot
break_sequence <- function(x, n = 6) {
  seq(min(x), max(x), length.out = n)
}

column_annotation_data <- tibble(
  pseudotime = progressions$cumpercentage,
  annotation = seu_oi$annotation[progressions$cell_id],
  nUMI = log10(seu_oi$nCount_RNA[progressions$cell_id])
)

cluster_colors2 <- cluster_colors

pseudotime_colors <- circlize::colorRamp2(break_sequence(column_annotation_data$pseudotime), RColorBrewer::brewer.pal(6, "Greys"))
umi_colors <- circlize::colorRamp2(break_sequence(column_annotation_data$nUMI), RColorBrewer::brewer.pal(6, "YlOrRd"))
column_annotation <- ComplexHeatmap::columnAnnotation(
  `Pseudotime` = anno_simple(
    column_annotation_data$pseudotime,
    pseudotime_colors
  ),
  `Celltype` = anno_simple(column_annotation_data$annotation, cluster_colors2),
  `# UMIs` = anno_simple(
    column_annotation_data$nUMI,
    umi_colors
  )
)
column_annotation_legends <- list(
  Legend(col_fun = pseudotime_colors, title = "Pseudotime", border = "#333333", direction = "horizontal"),
  Legend(col_fun = umi_colors, title = "# UMIs", border = "#333333", at = c(3, 4, 5), labels = c("1k", "10k", "100k"), direction = "horizontal"),
  Legend(at = names(cluster_colors2), title = "Celltype", border = "#333333", legend_gp = gpar(fill = cluster_colors2), ncol = 3)
)

heatmap_legend_param <- list(direction = "horizontal", border = "#333333", at = c(0, 1), labels = c("low", "high"), legend_width = unit(3, "cm"), title_position = "topcenter")
heatmap <- Heatmap(
  expression_matrix,
  cluster_rows = clust,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  top_annotation = column_annotation,
  col = rev(RColorBrewer::brewer.pal(10, "RdBu")[2:9]),
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = heatmap_legend_param,
  name = "Expression"
)
heatmap_list <- ComplexHeatmap::HeatmapList()
heatmap_list <- ComplexHeatmap::add_heatmap(heatmap, heatmap_list)

annotation_legends = c(
  column_annotation_legends
)

heatmap <- ComplexHeatmap::make_layout(
  heatmap_list,
  annotation_legend_list = annotation_legends,
  merge_legend = TRUE,
  heatmap_legend_side = "top",
  annotation_legend_side = "top"
)
save_plot(heatmap, experiment$results_file("heatmap"), width = 15)

