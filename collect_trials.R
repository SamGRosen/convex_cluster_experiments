library(tidyverse)

pattern = "CONVEX_CLUSTERING_43336712"
path = paste0("/cwork/sgr26/", pattern)

print(paste("Getting files for path", path))

paths_to_combine <- list.files(path, full.names = TRUE)

print(paste("Reading", length(paths_to_combine), "RDS files"))

combined_df <- paths_to_combine |>
  map_dfr(read_csv, .id = "trial")

csv_path <- paste0("~/convex_cluster_experiments/", pattern, ".csv")

print(paste("Saving csv to", csv_path))

write.csv(combined_df, csv_path, row.names=FALSE)
