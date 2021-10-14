


######

end_milestone_ids <- trajectory$milestone_network %>% filter(!(to %in% from)) %>% pull(to)

trajectory$milestone_network

# features_oi <- milestone_feature_importances %>%
#   filter(milestone_id %in% end_milestone_ids) %>% 
#   group_by(milestone_id) %>%
#   slice(1:40) %>%
#   pull(feature_id) %>%
#   as.character() %>%
#   unique()

markers <- seu_oi %>% FindMarkers(group.by = "annotation", ident.1 = "Spermatocyte Cyst Cell B", ident.2 = "Spermatocyte Cyst Cell A")

n_markers <- 50
features_oi <- c(
  markers %>% arrange(avg_logFC) %>% head(n_markers) %>% rownames(),
  markers %>% arrange(desc(avg_logFC)) %>% head(n_markers) %>% rownames()
)
row_split <- c(rep(c("Spermatocyte Cyst Cell A"), n_markers), rep("Spermatocyte Cyst Cell B", n_markers))

linearised <- linearise_cells(trajectory)
progressions <- linearised$progressions# %>% group_by(from, to) %>% sample_n(100, replace = TRUE) %>% arrange(cumpercentage)
expression_matrix <- t(as.matrix(seu_oi@assays$RNA@data[features_oi, progressions$cell_id]))

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
  pseudotime = calculate_pseudotime(trajectory)[progressions$cell_id],
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
  Legend(at = names(cluster_colors), title = "Celltype", border = "#333333", legend_gp = gpar(fill = cluster_colors2), ncol = 4)
)

column_split <- linearised$progressions$edge_id

heatmap_legend_param <- list(direction = "horizontal", border = "#333333", at = c(0, 1), labels = c("low", "high"), legend_width = unit(3, "cm"), title_position = "topcenter")
heatmap <- Heatmap(
  expression_matrix,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  top_annotation = column_annotation,
  col = rev(RColorBrewer::brewer.pal(10, "RdBu")[2:9]),
  row_names_gp = gpar(fontsize = 8),
  # heatmap_legend_param = heatmap_legend_param,
  name = "Expression",
  column_split = column_split,
  # row_split = row_split,
  column_title = NULL
)
heatmap_list <- ComplexHeatmap::HeatmapList()
heatmap_list <- ComplexHeatmap::add_heatmap(heatmap, heatmap_list)

annotation_legends = c(
  column_annotation_legends
)

heatmap <- ComplexHeatmap::make_layout(
  heatmap_list,
  # annotation_legend_list = annotation_legends,
  merge_legend = TRUE,
  heatmap_legend_side = "top",
  annotation_legend_side = "top"
)
save_plot(heatmap, experiment$results_file("heatmap"), width = 15)

