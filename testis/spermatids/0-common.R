clusters <- c(
  "early elongation stage spermatid",
  "early-mid elongation-stage spermatid",
  "mid-late elongation-stage spermatid"
)

cluster_colors <- setNames(c("#FEAF77FF", "#fee451ff", "#afe301ff"), clusters)

scale_cluster_color <- scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))
scale_cluster_fill <- scale_fill_manual(values = cluster_colors, breaks = names(cluster_colors))


rank_features <- function(expression_matrix, q = 0.9) {
  apply(
    expression_matrix,
    1,
    function(expression) {
      mean(which(expression > quantile(expression, q)))
    }
  )
}
