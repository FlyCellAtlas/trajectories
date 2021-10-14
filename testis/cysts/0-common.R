library(tidyverse)

cluster_colors <- c(
  "cyst stem cell" = "#d53e4f",
  "early cyst cells" = "#fc8d59",
  "Spermatocyte Cyst Cell B" = "#fdc116",
  "Spermatocyte Cyst Cell A" = "#3288bd",
  "cyst cell branch" = "#3f8f39ff"
)

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