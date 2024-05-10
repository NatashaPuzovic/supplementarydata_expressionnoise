# script to generate founding network topologies
#############  Barabasi-Alberts networks

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

# tolerable level of density distance from specified, for generation
tol_dens_deviation <- 0.005

# all densities
ALLnetwork_densities <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

# number of edges added each step 
#numEdgesAddedEachStep <- ceiling(num_nodes/20) # ~ density 0.05
#numEdgesAddedEachStep <- ceiling(num_nodes/10) # ~ density 0.1
#numEdgesAddedEachStep <- ceiling(num_nodes/5) # ~ density 0.2
#numEdgesAddedEachStep <- ceiling(num_nodes/3)  # ~ density 0.3
#numEdgesAddedEachStep <- ceiling(num_nodes/2) # = density 0.4
#numEdgesAddedEachStep <- num_nodes - 1 # = density 0.5
#numEdgesAddedEachStep <- num_nodes*2 # this will be density 0.5 anyway
ALLnumEdgesAddedEachStep <- c(ceiling(num_nodes/20),
                              ceiling(num_nodes/10),
                              ceiling(num_nodes/5),
                              ceiling(num_nodes/3),
                              ceiling(num_nodes/2),
                              num_nodes - 1)

# this is for density = 0.05
param_ind = 1

#for(param_ind in 1:length(ALLnumEdgesAddedEachStep)){
  
  numEdgesAddedEachStep <- ALLnumEdgesAddedEachStep[param_ind]
  network_density <- ALLnetwork_densities[param_ind]
  
  print(paste(numEdgesAddedEachStep))
  
  # make dir
  if(!dir.exists(paste0("20211201_BA_", num_nodes, "_node", network_density, "dens"))){
    dir.create(paste0("20211201_BA_", num_nodes, "_node", network_density, "dens"))
    }
  setwd(paste0("20211201_BA_", num_nodes, "_node", network_density, "dens"))
  
  for(netNum in 1:num_networks){
    
    # generate a random directed network with one component (no unconnected genes) 
    # power = the power of the preferential attachment, the default is one, ie. linear preferential attachment.
    g <- barabasi.game(n = num_nodes, power = 2, directed = TRUE, out.pref = T, m = numEdgesAddedEachStep)
    while(components(g)$no != 1){
      g <- barabasi.game(n = num_nodes, power = 2, directed = TRUE, out.pref = T, m = numEdgesAddedEachStep)
    }
    
    print(paste(graph.density(g)))
    
    # convert to adj matrix & add weights
    rmtrx <- as_adjacency_matrix(g, sparse = FALSE) # directions flipped outward (only for BA!)
    #rvals <- rbinom(n = length(rmtrx[rmtrx != 0]), size = 1, prob = 0.5)
    #rvals[rvals == 0] <- -0.5
    #rvals[rvals == 1] <- 0.5
    #rmtrx[rmtrx != 0] <- rvals 
    rmtrx[rmtrx != 0] <- runif(n = length(which(rmtrx != 0) == TRUE), min = -3, max = 3)
    
    g <- graph_from_adjacency_matrix(t(rmtrx), weighted = TRUE, diag = TRUE, mode = "directed")

    # write out
    write.table(rmtrx, 
                file = sprintf("rmtrxFoundingNet_%d.txt", netNum), 
                sep = " ", row.names = FALSE, col.names = FALSE)
  
  }
#  setwd("..")
#}

##############################################  end   ############################################