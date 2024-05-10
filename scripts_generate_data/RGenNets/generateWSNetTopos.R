# script to generate founding network topologies
#############  Watts-Strogatz networks

rm(list=ls())

library(igraph)
library(ggplot2)

setwd("../../sim_data/")
if(!dir.exists("generated_net_topos")){
  dir.create("generated_net_topos")
}
setwd("generated_net_topos")
  
################### params

# number of genes in the network
num_nodes <- 40

# number of network topologies to create
num_networks <- 1000

# names of nodes
node_names <- paste("G", 1:num_nodes, sep = "")

# all densities
ALLnetwork_densities <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
ALLnei_size <- (1:7)*2

# this is for density = 0.05
param_ind = 1

#for(param_ind in 1:length(ALLnei_size)){
  
  nei_size <- ALLnei_size[param_ind]
  network_density <- ALLnetwork_densities[param_ind]
  
  # make dir
  if(!dir.exists(paste0("20211201_WS_", num_nodes, "_node", network_density, "dens"))){
    dir.create(paste0("20211201_WS_", num_nodes, "_node", network_density, "dens"))
    }
  
  setwd(paste0("20211201_WS_", num_nodes, "_node", network_density, "dens"))
  
  for(netNum in 1:num_networks){
    
    # generate a random directed network with one component (no unconnected genes) 
    # power = the power of the preferential attachment, the default is one, ie. linear preferential attachment.
    
    #if(nei_size == 1){
    #  g <- watts.strogatz.game(dim = 1, size = num_nodes, nei = paste(nei_size), p = paste(nei_size/20))
    #  while(components(g)$no != 1){
    #    g <- watts.strogatz.game(dim = 1, size = num_nodes, nei = paste(nei_size), p = paste(nei_size/20))
    #  }
      #plot(g)
    #} else {
    #  g <- watts.strogatz.game(dim = 1, size = num_nodes, nei = paste(nei_size), p = paste(nei_size/100))
    #  while(components(g)$no != 1){
    #    g <- watts.strogatz.game(dim = 1, size = num_nodes, nei = paste(nei_size), p = paste(nei_size/100))
    #  }
    #plot(g)
    #}
    
    g <- watts.strogatz.game(dim = 1, size = num_nodes, nei = paste(nei_size), p = paste(nei_size/20))
    while(components(g)$no != 1){
      g <- watts.strogatz.game(dim = 1, size = num_nodes, nei = paste(nei_size), p = paste(nei_size/20))
    }
    
    #print(paste(nei_size))
    #mean_distance(g)
    graph.density(g)
 
    # convert to adj matrix & add weights
    rmtrx <- as_adjacency_matrix(g, sparse = FALSE) # directions flipped outward
    #rvals <- rbinom(n = length(rmtrx[rmtrx != 0]), size = 1, prob = 0.5)
    #rvals[rvals == 0] <- -0.5
    #rvals[rvals == 1] <- 0.5
    #rmtrx[rmtrx != 0] <- rvals 
    rmtrx[rmtrx != 0] <- runif(n = length(which(rmtrx != 0) == TRUE), min = -3, max = 3)
    rmtrx[upper.tri(rmtrx)] <- 0 # make it directed
    
    g <- graph_from_adjacency_matrix(t(rmtrx), weighted = TRUE, diag = TRUE, mode = "directed")
    #plot(g, layout=layout.davidson.harel)
    print(graph.density(g))
    # write out
    write.table(rmtrx, 
                file = sprintf("rmtrxFoundingNet_%d.txt", netNum), 
                sep = " ", row.names = FALSE, col.names = FALSE)
  
  }
#  setwd("..")
#}
##############################################  end   ############################################