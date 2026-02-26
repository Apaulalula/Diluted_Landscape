###################################
# date: 02-2026
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Objective: Analyze different landscape configurations
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
###################################

######## Load Packages ######## 
library(betalink)
library(bipartite)
library(iNEXT)
library(vegan)
library(ggplot2)
library(network)
library(sna)
library(igraph)
library(ade4)
library(intergraph)
library(tibble)
library(dplyr)
library(tidyverse)
library(plyr)
library(reshape2)
library(cowplot)
library(ggpubr)

######## Load data ########

load("./Results_land_new/initial_p0.01_lnet_0.1.RData")
load("./Results_land_new/initial_p0.01_lnet_0.3.RData")
load("./Results_land_new/initial_p0.01_lnet_0.5.RData")
load("./Results_land_new/initial_p0.01_lnet_0.7.RData")
load("./Results_land_new/initial_p0.25_lnet_0.1.RData")
load("./Results_land_new/initial_p0.25_lnet_0.3.RData")
load("./Results_land_new/initial_p0.25_lnet_0.5.RData")
load("./Results_land_new/initial_p0.25_lnet_0.7.RData")
load("./Results_land_new/initial_p0.5_lnet_0.1.RData")
load("./Results_land_new/initial_p0.5_lnet_0.3.RData")
load("./Results_land_new/initial_p0.5_lnet_0.5.RData")
load("./Results_land_new/initial_p0.5_lnet_0.7.RData")

load("./Results_land_new/new_p0.01_lnet_0.1.RData")
load("./Results_land_new/new_p0.01_lnet_0.3.RData")
load("./Results_land_new/new_p0.01_lnet_0.5.RData")
load("./Results_land_new/new_p0.01_lnet_0.7.RData")
load("./Results_land_new/new_p0.25_lnet_0.1.RData")
load("./Results_land_new/new_p0.25_lnet_0.3.RData")
load("./Results_land_new/new_p0.25_lnet_0.5.RData")
load("./Results_land_new/new_p0.25_lnet_0.7.RData")
load("./Results_land_new/new_p0.5_lnet_0.1.RData")
load("./Results_land_new/new_p0.5_lnet_0.3.RData")
load("./Results_land_new/new_p0.5_lnet_0.5.RData")
load("./Results_land_new/new_p0.5_lnet_0.7.RData")


######################## Topology Measures ###################################
Net_metric_new <- data.frame() 

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    for (l in c("initial", "new")) {
      
    # Include meta_web in the list of local webs
    
    result_matrices <- get(paste0(l,"_p", p, "_lnet_", f))
    
    ######## Connectance #########
    
    connect_web <- sapply(result_matrices, networklevel, index="connectance")
    
    connect_web <- as.data.frame(connect_web)
    connect_web <- connect_web %>%
      rownames_to_column() %>% 
      separate(rowname, into = c("Local", "Metric"), sep = "\\.") %>% 
      select(-Metric)
    
    
    ######## Make graph object ######## 
    local_graphs <- list()
    
    for (i in unique(names(result_matrices))) {
      
      net_graph = graph_from_biadjacency_matrix(result_matrices[[i]], weighted = NULL)
      local_graphs[[i]]<- net_graph
      
    }
    
    ########## Make connected graphs #########
    net_connect <- list()
    
    for (c in unique(names(local_graphs))) {
      
      net_1<-(local_graphs[[c]])
      
      # Find the components of the graph
      components <- components(net_1)
      
      # Extract component membership and sizes
      membership <- components$membership
      sizes <- components$csize
      
      # Print components
      print(sizes)
      
      if(is_connected(net_1) == FALSE){
        # Connect the components by adding an edge
        # Find representative nodes from each component
        component_ids <- unique(membership)
        
        # Connect each pair of components by adding an edge between them
        for (d in 1:(length(component_ids) - 1)) {
          # Find a node from the current component
          node_from <- which(membership == component_ids[d])[1]
          # Find a node from the next component
          node_to <- which(membership == component_ids[d + 1])[1]
          
          # Add an edge between these nodes
          net_1 <- add_edges(net_1, c(node_from, node_to))
        }
        
      }
      # Check if the graph is now connected
      print(is_connected(net_1))  # Should return TRUE
      
      net_connect[[c]]<-net_1
      
    }
    
    ######## Size ######## 
    network_sizes <- lapply(local_graphs, vcount)
    network_edge_counts <- lapply(local_graphs, ecount)
    
    head(network_sizes)
    network_sizes = unlist(network_sizes)
    
    network_sizes <- as.data.frame(network_sizes)
    network_sizes <- network_sizes %>% 
      rownames_to_column() %>% 
      dplyr::rename("Local"= rowname)
    
    ######## Modularity ######## 
    mod.groups <- sapply(net_connect, cluster_fast_greedy, weights = NULL)
    like.m.groups <- sapply(mod.groups, modularity)
    
    like.m.groups <- data.frame(Modularity = like.m.groups)
    like.m.groups <- like.m.groups %>% 
      rownames_to_column() %>% 
      dplyr::rename("Local"= rowname)
    
    ######## Nestedness ######## 
    
    nest.groups <- sapply(result_matrices, nestednodf)
    nest.groups= as_tibble(nest.groups)
    nest.g = nest.groups[3,c(1:length(result_matrices))]  %>% unnest(cols = c(1:length(result_matrices)))
    
    nest.g<- as.data.frame(t(nest.g))
    colnames(nest.g)= c("N.collums", "N.rows", "NODF")
    
    nest.g <- nest.g %>% 
      mutate("Local"= as.character(unique(names(local_graphs))))
    
    ######## Creating Data frame with all metrics for all networks ######## 
    
    groups.metric = full_join(network_sizes, nest.g, by = "Local")
    groups.metric = full_join(groups.metric, like.m.groups, by = "Local")
    groups.metric = full_join(groups.metric, connect_web, by = "Local")
    
    colnames(groups.metric) = c("Local","Size", "N.collums", "N.rows", "NODF", "Modularity", "Connectance")
    
    groups.metric= groups.metric %>% 
      select(Local, Size, NODF, Modularity, Connectance) %>% 
      mutate(p_value = p,
             forest_cover = f,
             land = l)
    
    Net_metric_new <- rbind(Net_metric_new, groups.metric)
    }
  }
}

Net_metric_new1 <- Net_metric_new %>% 
  separate(Local, c(NA, "patch_no")) %>% 
  mutate(patch_no = as.integer(patch_no)) %>% 
  inner_join(land_pv) %>% 
  select(patch_no, forest_cover, p_value, land, Land_use, Size, NODF, Modularity, Connectance) %>% 
  dplyr::rename(patch_id = patch_no)

str(Net_metric_new1)

############### Summarise metrics ###########################
new_net_mean <- Net_metric_new %>% 
  separate(Local, c(NA, "patch_no")) %>% 
  mutate(patch_no = as.integer(patch_no)) %>% 
  inner_join(land_pv) %>% 
  select(patch_no, forest_cover, p_value, land, Land_use, Size, NODF, Modularity, Connectance) %>% 
  dplyr::group_by(forest_cover, p_value, land, Land_use) %>% 
  dplyr::summarise(Mean_connectance = mean(Connectance, na.rm = TRUE),
                   sd_con = sd(Connectance),
                   Mean_size = mean(Size, na.rm = TRUE),
                   sd_size= sd(Size),
                   Mean_NODF = mean(NODF, na.rm = TRUE),
                   sd_NODF = sd(NODF),
                   Mean_mod = mean(Modularity, na.rm = TRUE),
                   sd_mod= sd(Modularity),
                   .groups = "drop") %>% 
  mutate(zc_con = Mean_connectance/sd_con, zc_size = Mean_size/sd_size, zc_NODF = Mean_NODF/sd_NODF, zc_mod = Mean_mod/sd_mod)

write.csv(new_net_mean, "./Results_land_new/new_net_mean.csv")

