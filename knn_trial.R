source("./helpers.R")

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)

print(args)

stopifnot(length(args) == 2)

JOB_ID <- args[1]
ARRAY_ID <- as.integer(args[2])

set.seed(ARRAY_ID)

means <- matrix(c( 1,  1,
                   -1,  1,
                   1, -1,
                   -1, -1),
                nrow = 4, ncol = 2, byrow = T)
means <- means * sqrt(2)

cov_1 <- matrix(c(1, 0,
                  0, 1), nrow = 2)

covs <- array(c(cov_1, cov_1, cov_1, cov_1), dim = c(2, 2, 4))

trial = generate_trial(4, 25, means, covs)

gammas <- seq(0, 100, 0.01)

phi <- 1

make_row_for_run <- function(run, method, extra_params = "") {
  eff_res_df <- as.data.frame(run$pseudoinverse_info$effective_resistance) |>
    mutate(from = row_number()) |>
    pivot_longer(!from, names_to = "to", values_to = "eff_res") |>
    mutate(to = as.numeric(sub("V", "", to))) |>
    mutate(from_label = run$trial$labels[from],
           to_label = run$trial$labels[to]) |>
    inner_join(
      as.data.frame(as.matrix(dist(run$trial$X))) |>
        mutate(from = row_number()) |>
        pivot_longer(!from, names_to = "to", values_to = "dist") |>
        mutate(to = as.numeric(sub("V", "", to))),
      by = c("from", "to")
    ) |>
    inner_join(
      as.data.frame(as.matrix(dist(run$trial$U))) |>
        mutate(from = row_number()) |>
        pivot_longer(!from, names_to = "to", values_to = "true_dist") |>
        mutate(to = as.numeric(sub("V", "", to))),
      by = c("from", "to")
    )

  if(ncol(run$pseudoinverse_info$F_dagger) != nrow(run$trial$X)) {
    F_dagger_E <- NA
    max_F_dagger_E_row <- NA
    sse_for_theorem_gamma <- NA
    ari_for_theorem_gamma <- NA
  } else {
    F_dagger_E <- run$pseudoinverse_info$F_dagger %*% (run$trial$X - run$trial$U)
    max_F_dagger_E_row <- max(sqrt(rowSums(F_dagger_E^2)))

    sse_for_theorem_gamma <- run$cluster_sse |>
      filter(lambda > max_F_dagger_E_row) |>
      slice_min(lambda) |>
      pull(sse)

    ari_for_theorem_gamma <- run$cluster_ari |>
      filter(max_F_dagger_E_row >= smallest_lambda,
             max_F_dagger_E_row <= largest_lambda) |>
      pull(ari)

    if(length(sse_for_theorem_gamma) == 0) {
      sse_for_theorem_gamma = NA
    }

    if(length(ari_for_theorem_gamma) == 0) {
      ari_for_theorem_gamma = NA
    }
  }

  incidence_rank <- rankMatrix(run$pseudoinverse_info$F_dagger)

  tibble_row(
    method = method,
    extra_params = extra_params,
    best_ari = run$cluster_ari |> slice_max(ari, with_ties = F) |> pull(ari),
    best_ari_num_cluster = run$cluster_ari |> slice_max(ari, with_ties = F) |> pull(num_cluster),
    best_ari_smallest_lambda = run$cluster_ari |> slice_max(ari, with_ties = F) |> pull(smallest_lambda),
    best_ari_largest_lambda = run$cluster_ari |> slice_max(ari, with_ties = F) |> pull(largest_lambda),
    best_sse = run$cluster_sse |> slice_min(sse, with_ties = F) |> pull(sse),
    best_sse_lambda = run$cluster_sse |> slice_min(sse, with_ties = F) |> pull(lambda),
    within_cluster_edges = sum(run$plot_dfs$edges$label == run$plot_dfs$edges$labelend),
    between_cluster_edges = sum(run$plot_dfs$edges$label != run$plot_dfs$edges$labelend),
    smallest_clusters_fit = (run$clusterpath$clusters_df |> slice_min(num_cluster) |> pull(num_cluster))[1],
    within_cluster_eff_resistance = sum(eff_res_df |> filter(from_label == to_label) |> pull(eff_res)),
    between_cluster_eff_resistance = sum(eff_res_df |> filter(from_label != to_label) |> pull(eff_res)),
    within_cluster_dist = sum(eff_res_df |> filter(from_label == to_label) |> pull(dist)),
    between_cluster_dist = sum(eff_res_df |> filter(from_label != to_label) |> pull(dist)),
    between_cluster_true_dist = sum(eff_res_df |> filter(from_label != to_label) |> pull(true_dist)),
    within_cluster_eff_dist_product = sum((eff_res_df |> filter(from_label == to_label) |> pull(eff_res)) *
                                            (eff_res_df |> filter(from_label == to_label) |> pull(dist))),
    between_cluster_eff_dist_product = sum((eff_res_df |> filter(from_label != to_label) |> pull(eff_res)) *
                                             (eff_res_df |> filter(from_label != to_label) |> pull(dist))),
    within_cluster_F_dagger_row_norm = sum(run$plot_dfs$edges |> filter(label == labelend) |> pull(F_dagger_row_norm)),
    between_cluster_F_dagger_row_norm = sum(run$plot_dfs$edges |> filter(label != labelend) |> pull(F_dagger_row_norm)),
    max_F_dagger_row_norm = max(run$plot_dfs$edges$F_dagger_row_norm),
    F_dagger_E_max = max(abs(F_dagger_E)),
    F_dagger_E_mean_abs = mean(abs(F_dagger_E)),
    incidence_rank = incidence_rank,
    max_F_dagger_E_row = max_F_dagger_E_row,
    sse_for_theorem_gamma = sse_for_theorem_gamma,
    ari_for_theorem_gamma = ari_for_theorem_gamma
  )
}

to_save <- tibble()

for(k in 1:(nrow(trial$X) - 1)) {
  knn_weights <- sparse_weights(trial$X, k, 0, connected = F, scale = F)
  knn_weights_run <- get_all_info(trial, knn_weights, gammas)
  to_save <- to_save |>
    bind_rows(make_row_for_run(knn_weights_run, "knn", paste0("k=", k)))

  knn_weights_kernel <- sparse_weights(trial$X, k, phi, connected = F, scale = F)
  knn_weights_kernel_run <- get_all_info(trial, knn_weights_kernel, gammas)
  to_save <- to_save |>
    bind_rows(make_row_for_run(knn_weights_kernel_run, "knn_weighted", paste0("k=", k)))
}

output_dir <- paste0(
  "~/CONVEX_CLUSTERING_KNN_",
  JOB_ID,
  "/"
)

if (!dir.exists(output_dir)) {
  print(paste("Creating dir", output_dir))
  dir.create(output_dir)
} else {
  print(paste("Saving to existing", output_dir))
}

save_success <- tryCatch(
  expr = {
    write.csv(
      to_save,
      file = paste0(output_dir, "CONVEX_CLUSTERING_KNN_", JOB_ID, "_", ARRAY_ID, ".csv"),
      row.names = F)
    TRUE
  },
  error = function(e) {
    print(e)
    FALSE
  }
)

if (!save_success) {
  q(save = "no",
    status = 11,
    runLast = FALSE)
}
