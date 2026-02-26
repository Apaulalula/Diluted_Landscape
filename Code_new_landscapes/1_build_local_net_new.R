###################################
# date: 02-2026
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Objective: Analyze different landscape configurations
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
###################################

# Load Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(bipartite)

# Combining data
host_info <- read.csv("./Results/host_info.csv") %>% column_to_rownames(var = "X")
prst_info <- read.csv("./Results/prst_info.csv") %>% column_to_rownames(var = "X")

# gather model data
Sensitivity_land <- data.frame(forest_cover=double(),
                                p_value=double(),
                                land= double(),
                                patch_no=integer(),
                                guild= vector(),
                                species = vector(),
                                z= double())


for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
      
      l <- read.csv(paste0("./Output_v2b/metanetwork_fc",f,"_p",p,"_er1_ec0.5_cr1_cc0.6_mr0.5_mc0.7_phi0.5_alpha0.2_eps5.csv"), quote="\"", comment.char="")
      l <- data.frame(forest_cover=f,
                      p_value=p,
                      land= "initial",
                      patch_no=as.vector(l$patch_no),
                      guild=as.vector(l$guild),
                      species = as.vector(l$species),
                      z= as.vector(l$z))
      
      Sensitivity_land <- rbind(Sensitivity_land, l)
  }
}

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    
    l <- read.csv(paste0("./Output_v2b_new/metanetwork_fc",f,"_p",p,"_er1_ec0.5_cr1_cc0.6_mr0.5_mc0.7_phi0.5_alpha0.2_eps5.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    land= "new",
                    patch_no=as.vector(l$patch_no),
                    guild=as.vector(l$guild),
                    species = as.vector(l$species),
                    z= as.vector(l$z))
    
    Sensitivity_land <- rbind(Sensitivity_land, l)
  }
}

glimpse(Sensitivity_land)

save(Sensitivity_land,file="./Results_new/Sensitivity_land.RData")

############### Species data ###############
host_df <- host_info %>%
  rownames_to_column(var= "species") %>% 
  mutate(guild = "resources")

prst_df <- prst_info %>% 
  mutate(guild = "consumers",
         ZoonoticStatus = case_when(
           ZoonoticStatus == 0 ~ "Non-zoonotic",
           ZoonoticStatus == 1 ~ "Zoonotic"))

str(Sensitivity_land)

prst_df <- prst_df %>%  select(-X)

land_p <- Sensitivity_land %>%
      filter(guild == "consumers") %>%
      full_join(land_pv) %>%
      left_join(prst_df) %>% 
      select(forest_cover,p_value,land, patch_no, Land_use, Name, species, z, ZoonoticStatus)

land_h <- Sensitivity_land %>%
  filter(guild == "resources") %>% 
  full_join(land_pv) %>%
  mutate(species = as.character(species)) %>% 
  inner_join(host_df) %>% 
  select(forest_cover,p_value,land,patch_no, Land_use, Name, z, pf)

unique(land_p$Name)

#Build Local networks
# Get the names of host and parasite species
Minc <- as.matrix(Minc)

host_species <- rownames(Minc)
parasite_species <- colnames(Minc)

for(p in c(0.01,0.25,0.5)){
  for(f in c(0.1,0.3,0.5,0.7)){
    for (l in c("initial", "new")) {
    
    result_hosts <-  land_h %>% 
      filter(p_value == p) %>% 
      filter(forest_cover == f) %>% 
      filter(land == l)
    
    result_parasites <-land_p %>% 
      filter(p_value == p) %>% 
      filter(forest_cover == f) %>% 
      filter(land == l)
    
    # Number of sites
    num_sites <- unique(result_parasites$patch_no)
    
    # List to store the interaction matrices for each site
    var_name <- paste0(l,"_p", p, "_lnet_", f)
    
    site_list <- vector("list", length(num_sites))
    
    names(site_list) <- paste0("p_", num_sites)
    
    # Loop over each site
    for (site in num_sites) {
      print(site)
      # # Get the presence data for this site
      host_present <- result_hosts[result_hosts$patch_no == site, ]
      parasite_present <- result_parasites[result_parasites$patch_no == site, ]
      
      current_incidence <- as.matrix(Minc[host_present$Name, parasite_present$Name])
      current_incidence <- as.matrix(current_incidence[,colSums(current_incidence)>0, drop = FALSE])
      current_incidence <- as.matrix(current_incidence[rowSums(current_incidence)>0, , drop = FALSE])
      
      # Store the site-specific interaction matrix in the list
      
      site_list[[paste0("p_", site)]] <- current_incidence
      
    }
    
    assign(var_name, site_list)
    
    }
  }
}


for(p in c(0.01,0.25,0.5)){
  for(f in c(0.1,0.3,0.5,0.7)){
    for (l in c("initial", "new")) {
    
    var_name <- paste0(l,"_p", p,"_lnet_", f)
    
    filtered_matrices <- lapply(get(paste0(l,"_p", p,"_lnet_", f)), function(mat) {
      if (nrow(mat) > 2 && ncol(mat) > 2) {
        return(mat)
      } else {
        return(NULL)
      }
    })
    
    # Remove NULL elements (i.e., matrices that were too small)
    filtered_matrices <- Filter(Negate(is.null), filtered_matrices)
    assign(var_name, filtered_matrices)
    }
  }
}

save(initial_p0.01_lnet_0.1,file="./Results_land_new/initial_p0.01_lnet_0.1.RData")
save(initial_p0.01_lnet_0.3,file="./Results_land_new/initial_p0.01_lnet_0.3.RData")
save(initial_p0.01_lnet_0.5,file="./Results_land_new/initial_p0.01_lnet_0.5.RData")
save(initial_p0.01_lnet_0.7,file="./Results_land_new/initial_p0.01_lnet_0.7.RData")

save(initial_p0.25_lnet_0.1,file="./Results_land_new/initial_p0.25_lnet_0.1.RData")
save(initial_p0.25_lnet_0.3,file="./Results_land_new/initial_p0.25_lnet_0.3.RData")
save(initial_p0.25_lnet_0.5,file="./Results_land_new/initial_p0.25_lnet_0.5.RData")
save(initial_p0.25_lnet_0.7,file="./Results_land_new/initial_p0.25_lnet_0.7.RData")

save(initial_p0.5_lnet_0.1,file="./Results_land_new/initial_p0.5_lnet_0.1.RData")
save(initial_p0.5_lnet_0.3,file="./Results_land_new/initial_p0.5_lnet_0.3.RData")
save(initial_p0.5_lnet_0.5,file="./Results_land_new/initial_p0.5_lnet_0.5.RData")
save(initial_p0.5_lnet_0.7,file="./Results_land_new/initial_p0.5_lnet_0.7.RData")


save(new_p0.01_lnet_0.1,file="./Results_land_new/new_p0.01_lnet_0.1.RData")
save(new_p0.01_lnet_0.3,file="./Results_land_new/new_p0.01_lnet_0.3.RData")
save(new_p0.01_lnet_0.5,file="./Results_land_new/new_p0.01_lnet_0.5.RData")
save(new_p0.01_lnet_0.7,file="./Results_land_new/new_p0.01_lnet_0.7.RData")

save(new_p0.25_lnet_0.1,file="./Results_land_new/new_p0.25_lnet_0.1.RData")
save(new_p0.25_lnet_0.3,file="./Results_land_new/new_p0.25_lnet_0.3.RData")
save(new_p0.25_lnet_0.5,file="./Results_land_new/new_p0.25_lnet_0.5.RData")
save(new_p0.25_lnet_0.7,file="./Results_land_new/new_p0.25_lnet_0.7.RData")

save(new_p0.5_lnet_0.1,file="./Results_land_new/new_p0.5_lnet_0.1.RData")
save(new_p0.5_lnet_0.3,file="./Results_land_new/new_p0.5_lnet_0.3.RData")
save(new_p0.5_lnet_0.5,file="./Results_land_new/new_p0.5_lnet_0.5.RData")
save(new_p0.5_lnet_0.7,file="./Results_land_new/new_p0.5_lnet_0.7.RData")

# Initialize an empty list to store interaction data frames
interaction_list <- list()
idx <- 1  # Index for tracking list position

# Loop through combinations of alpha, p_value, and forest_cover
  for(f in c(0.1, 0.3, 0.5, 0.7)) {
    for(p in c(0.01, 0.25, 0.5)) {
      for (l in c("initial", "new")) {
      
      # Get the list of patches
      local_net <- get(paste0(l, "_p", p, "_lnet_", f))
      unique_patch <- names(local_net)
      
      # Loop through each patch matrix
      for(patch in unique_patch) {
        
        # Convert the interaction matrix to a data frame and add metadata
        df <- as.data.frame(as.table(local_net[[patch]])) %>% 
          mutate(patch_id = sub(".*_(.*)", "\\1", patch), 
                 forest_cover = f, 
                 p_value = p, 
                 land = l) %>% 
          dplyr::rename(host_species = Var1, 
                        parasite_species = Var2, 
                        value = Freq) %>% 
          filter(value == 1)
        
        # Store the data frame in the list
        interaction_list[[idx]] <- df
        idx <- idx + 1
      }
    }
  }
}


# Combine all data frames into one
interaction_new_df <- dplyr::bind_rows(interaction_list)

# Save data frames

write.csv(interaction_new_df, "./Results_land_new/interaction_new_df.csv")
