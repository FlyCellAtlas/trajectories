library(tidyverse)

cluster_colors <- c(
  "cyst stem cell" = "#d53e4f",
  "early cyst cell 1" = "#fc8d59",
  "early cyst cell 2" = "#fdc116",
  "spermatocyte cyst cell branch a" = "#3288bd",
  "spermatocyte cyst cell branch b" = "#b3d594ff",
  "cyst cell branch b" = "#3f8f39ff",
  "late cyst cell branch a" = "#3288bd",
  "late cyst cell branch b" = "#81bbdd",
  "terminal epithelial cell of testis" = "black"
)

cluster_colors <- c(
  "cyst stem cell" = "#d53e4f",
  "early cyst cell 1" = "#fc8d59",
  "early cyst cell 2" = "#fdc116",
  "spermatocyte cyst cell branch a" = "#3288bd",
  "spermatocyte cyst cell branch b" = "#b3d594ff",
  "cyst cell branch b" = "#3f8f39ff",
  "late cyst cell branch a" = "#3288bd",
  "late cyst cell branch b" = "#81bbdd",
  "terminal epithelial cell of testis" = "black"
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