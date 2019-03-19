#!/usr/local/bin/Rscript

task <- dyncli::main()

library(jsonlite)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Hotfix for drop = FALSE problem in URD
URD:::floodPseudotimeCalc %>%
  deparse() %>%
  gsub("cells.visited], 1, combine.probs)", "cells.visited, drop = FALSE], 1, combine.probs)", ., fixed = TRUE) %>%
  parse(text = .) %>%
  eval(envir = environment(URD:::floodPseudotimeCalc)) %>%
  utils::assignInNamespace("floodPseudotimeCalc", ., ns = "URD")

library(URD)

#   ____________________________________________________________________________
#   Load data                                                               ####

count <- task$counts
parameters <- task$parameters
start_id <- task$priors$start_id

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
if (parameters$sigma.use <= 0) {
  parameters$sigma.use <- NULL
}
if (parameters$knn <= 0) {
  parameters$knn <- destiny:::find_dm_k(nrow(count))
  if (parameters$knn >= nrow(count)) {
    parameters$knn <- round(log10(nrow(count)) * 10)
  }
}

# just load the data, filtering has already been done
urd <- createURD(count.data = t(count), min.cells = 0, min.counts = 0, min.genes = 0)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# Calculate Diffusion Map
urd <- calcDM(
  urd,
  knn = parameters$knn,
  sigma.use = parameters$sigma.use,
  distance = parameters$distance
)

# Then we run 'flood' simulations
urd.floods <- floodPseudotime(
  urd,
  root.cells = start_id,
  n = parameters$n_floods,
  minimum.cells.flooded = 0,
  verbose = FALSE
)

# The we process the simulations into a pseudotime
urd <- floodPseudotimeProcess(
  urd,
  urd.floods,
  floods.name = "pseudotime",
  stability.div = parameters$stability.div
)

# Calculate PCA and tSNE
urd <- calcPCA(
  urd,
  mp.factor = parameters$mp.factor
)
urd <- calcTsne(
  urd,
  perplexity = parameters$perplexity,
  theta = parameters$theta,
  max_iter = parameters$max_iter
)

# Calculate graph clustering of these cells
urd <- graphClustering(
  urd,
  num.nn = parameters$num.nn,
  do.jaccard = parameters$do.jaccard,
  method = "Louvain"
)
cluster_name <- paste0("Louvain-", parameters$num.nn)

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
urd.ptlogistic <- pseudotimeDetermineLogistic(
  urd,
  "pseudotime",
  optimal.cells.forward = parameters$optimal.cells.forward,
  max.cells.back = parameters$max.cells.back,
  do.plot = FALSE
)

# Bias the transition matrix acording to pseudotime
urd.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(
  urd,
  "pseudotime",
  logistic.parameters = urd.ptlogistic
))

# Simulate the biased random walks from each tip
urd.walks <- simulateRandomWalksFromTips(
  urd,
  tip.group.id = cluster_name,
  root.cells = start_id,
  transition.matrix = urd.biased.tm,
  n.per.tip = parameters$n.per.tip,
  root.visits = parameters$root.visits,
  max.steps = 5000,
  verbose = FALSE
)

# Process the biased random walks into visitation frequencies
urd <- processRandomWalksFromTips(
  urd,
  urd.walks,
  n.subsample = parameters$n.subsample,
  verbose = FALSE
)

# Load the cells used for each tip into the URD object
urd.tree <- loadTipCells(urd, cluster_name)

tips.use <- unique(urd.tree@group.ids[,cluster_name])

# Build the tree
urd.tree <- buildTree(
  urd.tree,
  pseudotime = "pseudotime",
  tips.use = tips.use,
  divergence.method = parameters$divergence.method,
  cells.per.pseudotime.bin = parameters$cells.per.pseudotime.bin,
  bins.per.pseudotime.window = parameters$bins.per.pseudotime.window,
  p.thresh = parameters$p.thresh,
  dendro.cell.jitter = 0,
  dendro.cell.dist.to.tree = 0,
  save.all.breakpoint.info = TRUE
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Save trajectory                                                         ####

tree_layout <- urd.tree@tree$tree.layout
subs <- tree_layout %>% filter(do.mean)
tree_layout <- tree_layout %>%
  filter(!do.mean) %>%
  left_join(subs %>% select(node.1 = node.2, node.new = node.1), by = "node.1") %>%
  mutate(node.1 = ifelse(is.na(node.new), node.1, node.new)) %>%
  select(from = node.1, to = node.2, x1, x2, y1, y2)

cell_layout <- urd.tree@tree$cell.layout

comb <- crossing(cell_layout, tree_layout) %>%
  filter(x1 <= x & x <= x2 & y1 <= y & y <= y2) %>%
  group_by(cell) %>%
  slice(1) %>%
  ungroup()

progressions <- comb %>%
  mutate(percentage = (y - y1) / (y2 - y1)) %>%
  select(cell_id = cell, from, to, percentage)

# collect milestone network
milestone_network <- tree_layout %>%
  mutate(
    length = abs(y1 - y2),
    directed = FALSE
  ) %>%
  select(from, to, length, directed)


output <-
  dynwrap::wrap_data(cell_ids = urd.tree@tree$cell.layout$cell) %>%
  dynwrap::add_trajectory(
    milestone_network = milestone_network,
    progressions = progressions
  )  %>%
  dynwrap::add_timings(timings = checkpoints)

dyncli::write_output(output, task$output)
