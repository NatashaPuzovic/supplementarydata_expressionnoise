# script to get network metrics for each gene

rm(list = ls())

#setwd("/home/npuzovic/expressionnoise/results/clustRun0506/net_1")
#source("functions_for_all_scripts.R")

if(file.exists("netMetrics.txt")){
  stop("Net metrics file already present.")
}

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(statnet))
suppressPackageStartupMessages(library(intergraph))

# read established network topology
rmtrx <- as.matrix(read.table("rmtrxEstNet.txt", sep = " ", header = FALSE))
rmtrx <- rmtrx[, -dim(rmtrx)[2]]

# number of nodes in the network
num_nodes <- dim(rmtrx)[1]

# name each node 
node_names = paste("G", 1:num_nodes, sep = "")


# construct igraph object from a regulatory matrix
makeiGraphObjectFromRmtrx <- function(rmtrx) {
  
  # make adjacency matrix
  adjmtrx <- t(rmtrx)
  
  # construct igraph object
  network_object <- igraph::graph_from_adjacency_matrix(adjmtrx, mode = c("directed"), weighted = TRUE,
                                                diag = TRUE)
  # name nodes
  V(network_object)$name <- node_names
  
  return(network_object)
}


# make network object
net_igr_raw <- makeiGraphObjectFromRmtrx(rmtrx = rmtrx)
net_igr_abs <- makeiGraphObjectFromRmtrx(rmtrx = abs(rmtrx))

# convert igraph to statnet object
net_stat_raw <- intergraph::asNetwork(net_igr_raw)



############################################################

# make df with all metrics
net_metrics <- data.frame(
  
  # gene-specific metrics
  
  # with negative and positive weights
  "node_names" = node_names,
  "k_all_inclps" = igraph::centr_degree(net_igr_raw, mode = "all", loops = TRUE)$res,
  "k_in_inclps" = igraph::centr_degree(net_igr_raw, mode = "in", loops = TRUE)$res,
  "k_out_inclps" = igraph::centr_degree(net_igr_raw, mode = "out", loops = TRUE)$res,
  "k_all_exclps" = igraph::centr_degree(net_igr_raw, mode = "all", loops = FALSE)$res,
  "k_in_exclps" = igraph::centr_degree(net_igr_raw, mode = "in", loops = FALSE)$res,
  "k_out_exclps" = igraph::centr_degree(net_igr_raw, mode = "out", loops = FALSE)$res,
  "clo_all" = igraph::centr_clo(net_igr_raw, mode = "all")$res, # only for undirected
  "betw" = igraph::centr_betw(net_igr_raw, directed = TRUE)$res,
  #"eigen_dir" = igraph::centr_eigen(net_igr_raw, directed = TRUE)["vector"],
  "eigen_undir" = igraph::centr_eigen(net_igr_raw, directed = FALSE)["vector"],
  "str_all_inclps" = igraph::strength(net_igr_raw, mode = "all", loops = TRUE),
  "str_in_inclps" = igraph::strength(net_igr_raw, mode = "in", loops = TRUE),
  "str_out_inclps" = igraph::strength(net_igr_raw, mode = "out", loops = TRUE),
  "str_all_exclps" = igraph::strength(net_igr_raw, mode = "all", loops = FALSE),
  "str_in_exclps" = igraph::strength(net_igr_raw, mode = "in", loops = FALSE),
  "str_out_exclps" = igraph::strength(net_igr_raw, mode = "out", loops = FALSE),
  "hub_score" = igraph::hub_score(net_igr_raw)$vector,
  "auth_incwght" = igraph::authority_score(net_igr_raw)$vector,
  "auth_excwght" = igraph::authority_score(net_igr_raw, weights = NA)$vector,
  
  # absolute weights
  "sum_abs_inlinks" = rowSums(abs(rmtrx)),
  "sum_abs_outlinks" = colSums(abs(rmtrx)),
  "absstr_all_inclps" = igraph::strength(net_igr_abs, mode = "all", loops = TRUE),
  "absstr_in_inclps" = igraph::strength(net_igr_abs, mode = "in", loops = TRUE),
  "absstr_out_inclps" = igraph::strength(net_igr_abs, mode = "out", loops = TRUE),
  "absstr_all_exclps" = igraph::strength(net_igr_abs, mode = "all", loops = FALSE),
  "absstr_in_exclps" = igraph::strength(net_igr_abs, mode = "in", loops = FALSE),
  "absstr_out_exclps" = igraph::strength(net_igr_abs, mode = "out", loops = FALSE),
  
  # additional metrics from statnet
  "flow" = sna::flowbet(net_stat_raw),
  "load" = sna::loadcent(net_stat_raw),
  "info" = sna::infocent(net_stat_raw),
  "stress" = sna::stresscent(net_stat_raw),
  
  
  # whole network metrics
  "diam" = igraph::diameter(net_igr_raw, directed = TRUE, weights = NA), 
  "meandst" = igraph::mean_distance(net_igr_raw, directed = TRUE),
  "assort" = igraph::assortativity_degree(net_igr_raw, directed = TRUE),
  
  "cntr_degr_all" = igraph::centr_degree(net_igr_raw, mode = "all", loops = TRUE)$centralization,
  "tmax_degr_all" = igraph::centr_degree(net_igr_raw, mode = "all", loops = TRUE)$theoretical_max,
  "cntr_indegr" = igraph::centr_degree(net_igr_raw, mode = "in", loops = TRUE)$centralization,
  "tmax_indegr" = igraph::centr_degree(net_igr_raw, mode = "in", loops = TRUE)$theoretical_max,
  "cntr_outdegr" = igraph::centr_degree(net_igr_raw, mode = "out", loops = TRUE)$centralization,
  "tmax_outdegr" = igraph::centr_degree(net_igr_raw, mode = "out", loops = TRUE)$theoretical_max,
  "cntr_clo_all" = igraph::centr_clo(net_igr_raw, mode = "all")$centralization,
  "tmax_clo_all" = igraph::centr_clo(net_igr_raw, mode = "all")$theoretical_max,
  "cntr_betw" = igraph::centr_betw(net_igr_raw, directed = TRUE)$centralization,
  "tmax_betw" = igraph::centr_betw(net_igr_raw, directed = TRUE)$theoretical_max,
  "ave_k_all_inclps" = mean(igraph::centr_degree(net_igr_raw, mode = "all", loops = TRUE)$res),
  "ave_k_in_inclps" = mean(igraph::centr_degree(net_igr_raw, mode = "in", loops = TRUE)$res),
  "ave_k_out_inclps" = mean(igraph::centr_degree(net_igr_raw, mode = "out", loops = TRUE)$res),
  "ave_k_all_exclps" = mean(igraph::centr_degree(net_igr_raw, mode = "all", loops = FALSE)$res),
  "ave_k_in_exclps" = mean(igraph::centr_degree(net_igr_raw, mode = "in", loops = FALSE)$res),
  "ave_k_out_exclps" = mean(igraph::centr_degree(net_igr_raw, mode = "out", loops = FALSE)$res)
)

# ADD AND CHECK THIS
#colnames(evolnet)[c(22,23)] <- c("eigen_dir", "eigen_undir")

# output all stats
write.table(net_metrics, 
            file = "netMetrics.txt", 
            sep = " ", row.names = TRUE, col.names = TRUE)


############################################################        END         ########################################################