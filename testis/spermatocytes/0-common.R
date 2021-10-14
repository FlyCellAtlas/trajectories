library(tidyverse)
clusters <- c(
  "spermatogonium",
  "mid-late proliferating spermatogonia",
  "spermatogonium-spermatocyte transition",
  "spermatocyte",
  "late primary spermatocytes",
  "spermatid"
)

cluster_colors <- setNames(viridis::magma(length(clusters) + 1)[1:length(clusters)], clusters)

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
