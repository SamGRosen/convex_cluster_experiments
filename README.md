# Affinity Graph Connectivity in Convex Clustering Experiments

R code to reproduce figures and run all trials through slurm.

## Available Code

* `affinity_graph_illustration.qmd`: Produce figures of affinity graphs
* `affinity_solution_path.qmd`: Produce figures of affinity graphs and solution paths
* `barbell.qmd`: Draw an example ideal affinity graph
* `experiment_figures.qmd`: Produce figures summarizing all trial results
* `helpers.R`: A variety of documented functions used throughout this repository
* `knn_trial.R`: Generate a single dataset and run convex clustering for all k-nearest neighbor graphs
* `run_knn_trials.sh`: Run `knn_trial.R` through slurm
* `run_trials.sh`: Run `trial.R` through slurm
* `trial.R`: Generate a single dataset and run convex clustering over many randomly sampled affinity graphs


## Required packages

1. `tidyverse`
2. `CCMMR` 
3. `MASS`
4. `Matrix`
5. `aricode`
6. `gtools`
7. `ggrepel`
8. `scales`
9. `ggthemes`

## Reproducing results

1. Run Slurm
```
sbatch run_knn_trials.sh
sbatch run_trials.sh
```

2. Run `collect_trials.R` to consolidate the CSV files generated.
3. Run `experiment_figures.qmd` notebook to produce figures.
