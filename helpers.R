library(dplyr)
library(CCMMR)
library(MASS)
library(Matrix)
library(aricode)
library(gtools)
library(tidyr)
library(tibble)


#' Generate data from a Gaussian mixture model with given means and covariances.
#'
#' @param num_clusters Number of clusters
#' @param per_cluster Number of points to generate per cluster
#' @param means Matrix with cluster means as rows
#' @param covs Array with cluster covariance in last index
#'
#' @returns list with `X` as generated data, `U` as true centers, and `labels` as labels
generate_trial <- function(num_clusters, per_cluster, means, covs) {
  X <- matrix(NA, nrow = 0, ncol = ncol(means))
  for(cluster in 1:num_clusters) {
    X <- X |> rbind(mvrnorm(per_cluster, means[cluster,], covs[,,cluster]))
  }

  U <- matrix(NA, nrow = 0, ncol = ncol(X))
  for(cluster in 1:num_clusters) {
    U <- U |>
      rbind(matrix(means[cluster,], nrow = per_cluster, ncol = 2, byrow = T))
  }

  labels <- rep(1:num_clusters, each = per_cluster)
  list(
    X = X,
    U = U,
    labels = labels
  )
}


#' Generate affinities with all points in the same cluster having positive affinites
#' and a single affinity per pair of clusters
#'
#' @param per_cluster Size of each cluster
#' @param clusters Number of clusters
#'
#' @returns `sparseweights` object for input to `CCMMR`
get_oracle_graph <- function(per_cluster, clusters) {
  keys <- gtools::permutations(per_cluster, 2)
  for(cluster in 1:(clusters-1)) {
    keys <- rbind(keys, gtools::permutations(per_cluster, 2) + per_cluster * (cluster))
  }

  within <- rep(1, nrow(keys))
  cluster_pairs <- gtools::permutations(clusters, 2)
  extra_keys <- matrix(-1, nrow = nrow(cluster_pairs), ncol = 2)

  for(cluster_pair in 1:nrow(cluster_pairs)) {
    c1 <- cluster_pairs[cluster_pair, 1]
    c2 <- cluster_pairs[cluster_pair, 2]
    extra_keys[cluster_pair, 1] = 1 + (c1 - 1) * per_cluster
    extra_keys[cluster_pair, 2] = 1 + (c2 - 1) * per_cluster
  }


  keys <- rbind(keys, extra_keys)
  within <- c(within, rep(0, nrow(extra_keys)))

  edge_direction <- keys[,1] > keys[,2]
  keys <- keys[edge_direction, ]
  within <- within[edge_direction]

  values <- rep(1, nrow(keys))

  structure(
    list(values = values, keys = keys, within = within),
    class = "sparseweights"
  )
}


#' Generate random affinities with specific edge connecitivty within/between cluster
#'
#' @param per_cluster Number of points per cluster
#' @param clusters Number of clusters
#' @param within_density Percent of edges within cluster to generate
#' @param between_density Percent of edges between clusters to generate
#'
#' @returns `sparseweights` object for input to `CCMMR`
get_sparse_oracle_graph <- function(per_cluster, clusters,
                                    within_density = 1, between_density = 1) {
  keys <- matrix(NA, nrow = 0, ncol = 2)
  within <- c()
  for(cluster in 1:clusters) {
    cluster_keys <- 1:per_cluster + per_cluster * (cluster - 1)
    for(other_cluster in 1:clusters) {

      other_keys <- 1:per_cluster + per_cluster * (other_cluster - 1)

      combos <- expand.grid(cluster_keys, other_keys)

      if(cluster != other_cluster) {
        new_edges <- combos[runif(nrow(combos)) < between_density,]
      } else {
        new_edges <- combos[runif(nrow(combos)) < within_density,]
      }
      keys <- rbind(keys, new_edges)
      within <- c(within, rep(cluster == other_cluster, nrow(new_edges)))
    }
  }

  edge_direction <- keys[,1] > keys[,2]
  keys <- keys[edge_direction, ]
  within <- within[edge_direction]
  values <- rep(1, nrow(keys))

  structure(
    list(values = values, keys = keys, within = within),
    class = "sparseweights"
  )
}


#' Given data, generate Gaussian kernel affinities
#'
#' @param X Input data, with samples as rows
#' @param phi Kernel hyperparameter: exp(-phi * \|X_i - X_j\|_2^2)
#'
#' @returns `sparseweights` object for input to `CCMMR`
get_fully_connected_weighted_graph <- function(X, phi) {
  d <- dist(X)
  d2 <- as.vector(d)^2
  values <- exp(-d2 * phi)
  keys <- t(combn(nrow(X), 2))

  structure(
    list(values = values, keys = keys[, c(2, 1)]),
    class = "sparseweights"
  )
}


#' Generate dataframes useful for plotting affinity graphs
#'
#' @param X input data
#' @param labels labels for each data point
#' @param sparseweights `sparseweights` object used by CCMMR
#'
#' @returns list with `points` and `edges` keys for plotting
get_plot_dfs <- function(X, labels, sparseweights) {
  points_to_plot <- tibble(
    x = X[, 1],
    y = X[, 2],
    id = 1:nrow(X),
    label = labels
  )

  k1 <- sparseweights$keys[, 1]
  k2 <- sparseweights$keys[, 2]

  diffs <- X[k1, , drop = FALSE] - X[k2, , drop = FALSE]
  distances <- sqrt(rowSums(diffs * diffs))

  edges_as_df <- tibble(
    from = sparseweights$keys[, 1],
    to = sparseweights$keys[, 2],
    value = sparseweights$values,
    dist = distances
  )

  segments <- edges_as_df |>
    left_join(points_to_plot, by = c("from" = "id")) |>
    left_join(points_to_plot, by = c("to" = "id"), suffix = c("", "end")) |>
    mutate(
      `Within-component Edge` = case_when(
        label == labelend ~ "Yes",
        label != labelend ~ "No"
      )
    )

  list(
    points = points_to_plot,
    edges = segments
  )
}


#' Create a sparse incidence matrix from an edge dataset
#'
#' @param edge_dataset `edges` dataframe from `get_plot_dfs`
#'
#' @returns `sparseMatrix`
make_incidence_matrix <- function(edge_dataset) {
  graph_edges <- edge_dataset |> filter(
    from > to
  )

  num_edges <- nrow(graph_edges)
  sparseMatrix(
    i = c(graph_edges |> pull(from), graph_edges |> pull(to)),
    j = c(1:num_edges, 1:num_edges),
    x = c(-sqrt(graph_edges |> pull(value)),
          sqrt(graph_edges |> pull(value)))
  )
}


#' From an edge dataset, generate affinity graph info from the pseudoinverse of
#' the incidence matrix
#'
#' @param edge_dataset `edges` dataframe from `get_plot_dfs`
#'
#' @returns list with pseudoinverse as `F_dagger` key,
#'  `F_dagger_info` for diagnostics, and
#'  `effective_resistance` for effective resistance of each edge
get_pseudoinverse_info <- function(edge_dataset) {
  graph_edges <- edge_dataset |> filter(
    from > to
  )
  incidence_matrix <- make_incidence_matrix(edge_dataset)

  laplacian_pseudoinv <- ginv(as.matrix(crossprod(t(incidence_matrix))))
  L_dagger_diag <- base::diag(laplacian_pseudoinv)
  ones <- rep(1, nrow(laplacian_pseudoinv))

  effective_resistance <- outer(L_dagger_diag, ones, "*") +
                          outer(ones, L_dagger_diag, "*") -
                          2 * laplacian_pseudoinv

  F_dagger <- t(incidence_matrix) %*% laplacian_pseudoinv

  F_dagger_info <- tibble(
    from = c(graph_edges |> pull(from)),
    to = c(graph_edges |> pull(to)),
    F_dagger_row_norm = sqrt(rowSums(F_dagger ^ 2)),
    effective_resistance = diag(F_dagger %*% incidence_matrix)
  )

  list(
    F_dagger_info = F_dagger_info,
    effective_resistance = effective_resistance,
    F_dagger = F_dagger
  )
}


#' Given input data, affinities and hyperparameters, calculate solution path with
#' clustering info.
#'
#' @param X Input data
#' @param weights input weights
#' @param gammas grid of hyperparameters
#'
#' @returns list with `clusters_df` and `solutions_df`
get_solutions_info <- function(X, weights, gammas) {
  clusterpath = convex_clusterpath(
    X, weights, gammas,
    scale = F,
    center = F,
    save_clusterpath = TRUE,
    save_losses = TRUE
  )

  solutions_df <- tibble(
    id = rep(1:clusterpath$n, length(clusterpath$lambdas)),
    x = clusterpath$coordinates[, 1],
    y = clusterpath$coordinates[, 2],
    lambda = rep(clusterpath$lambdas, each = clusterpath$n)
  )

  clusters_df <- tibble(
    id = numeric(),
    num_cluster = numeric(),
    est_label = numeric(),
    smallest_lambda = numeric()
  )

  clusters_df <- clusterpath$info |>
    group_by(clusters) |>
    reframe(
      id = 1:clusterpath$n,
      est_label = clusters(clusterpath, first(clusters)),
      smallest_lambda = min(lambda),
      largest_lambda = max(lambda)
    ) |>
    rename(num_cluster = clusters)

  list(
    clusters_df = clusters_df,
    solutions_df = solutions_df
  )
}


#' Find which hyperparameter groups two data points together first
#'
#' @param cluster_df from `get_solutions_info`
#'
#' @returns dataframe of smallest hyperparameter and pairs of points
get_first_connections <- function(cluster_df) {
  cross_join(cluster_df, cluster_df) |>
    filter(id.x != id.y, est_label.x == est_label.y, num_cluster.x == num_cluster.y) |>
    group_by(id.x, id.y) |>
    summarise(smallest_lambda = min(smallest_lambda.x),
              most_clusters = max(num_cluster.x),
              .groups = "drop")
}


#' Given a trial and `sparseweights`, calculate the solution path and generate
#' all info needed for experiments
#'
#' @param trial output of `generate_trial`
#' @param weights `sparseweights` object for `CCMMR`
#' @param gammas hyperparameter grid
#'
#' @returns large list with all info
get_all_info <- function(trial, weights, gammas) {
  plot_dfs <- get_plot_dfs(trial$X, trial$labels, weights)
  pseudoinverse_info <- get_pseudoinverse_info(plot_dfs$edges)

  clusterpath <- get_solutions_info(trial$X, weights, gammas)
  cluster_ari <- clusterpath$clusters_df |>
    group_by(smallest_lambda, largest_lambda, num_cluster) |>
    summarise(ari = aricode::ARI(est_label, trial$labels), .groups = "drop")

  cluster_sse <- clusterpath$solutions_df |>
    group_by(lambda) |>
    summarise(sse = sum((x - trial$U[, 1])^2 +
                          (y - trial$U[, 2])^2), .groups = "drop")

  plot_dfs$edges <- plot_dfs$edges |>
    inner_join(pseudoinverse_info$F_dagger_info, by = c("from", "to"))
  list(
    plot_dfs = plot_dfs,
    trial = trial,
    weights = weights,
    gammas = gammas,
    clusterpath = clusterpath,
    cluster_ari = cluster_ari,
    cluster_sse = cluster_sse,
    pseudoinverse_info = pseudoinverse_info
  )
}
