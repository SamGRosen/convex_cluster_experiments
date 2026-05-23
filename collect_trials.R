library(tidyverse)

pattern = "CONVEX_CLUSTERING_KNN_43363704"
path = paste0("~/", pattern)

print(paste("Getting files for path", path))

paths_to_combine <- list.files(path, full.names = TRUE)

print(paste("Reading", length(paths_to_combine), "RDS files"))

combined_df <- map_dfr(
  paths_to_combine,
  ~ read_csv(.x, show_col_types = FALSE, progress = FALSE),
  .progress = TRUE
)

csv_path <- paste0("~/convex_cluster_experiments/", pattern, ".csv")

print(paste("Saving csv to", csv_path))

write.csv(combined_df, csv_path, row.names=FALSE)
