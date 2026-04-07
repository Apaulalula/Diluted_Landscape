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
host_info <- read.csv("./Results/host_info.csv") %>% 
  column_to_rownames(var = "X")
prst_info <- read.csv("./Results/prst_info.csv") %>% 
  column_to_rownames(var = "X")

Minc <- read.csv("./Data/M_inc.csv") %>% 
  column_to_rownames(var = "X") %>% 
  dplyr::rename('juquitiba-like_virus_strain_on576' = 'juquitiba.like_virus_strain_on576')

# gather model data
Sensitivity_trait <- data.frame(forest_cover=double(),
                                p_value=double(),
                                patch_no=integer(),
                                guild= vector(),
                                species = vector(),
                                z= double())


for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    
    l <- read.csv(paste0("./Output_v2b_new/ICT_metanetwork_fc",f,"_p",p,"_er1_ec0.5_cr1_cc0.6_mr0.5_mc0.7_phi0.5_alpha0.2_eps5.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=as.vector(l$patch_no),
                    guild=as.vector(l$guild),
                    species = as.vector(l$species),
                    z= as.vector(l$z))
    
    Sensitivity_trait <- rbind(Sensitivity_trait, l)
  }
}

glimpse(Sensitivity_trait)

save(Sensitivity_trait,file="./Results_new/Sensitivity_trait.RData")

############### Species data ###############
host_df <- host_info %>%
  rownames_to_column(var= "species") %>% 
  mutate(guild = "resources") %>% 
  glimpse()

prst_df <- prst_info %>% 
  mutate(guild = "consumers",
         ZoonoticStatus = case_when(
           ZoonoticStatus == 0 ~ "Non-zoonotic",
           ZoonoticStatus == 1 ~ "Zoonotic")) %>% 
  glimpse()


land_p <- Sensitivity_trait %>%
      filter(guild == "consumers") %>%
      full_join(land_ict) %>%
      left_join(prst_df) %>% 
      select(forest_cover,p_value, patch_no, Land_use, Name, species, z)

land_h <- Sensitivity_trait %>%
  filter(guild == "resources") %>% 
  full_join(land_ict) %>%
  mutate(species = as.character(species)) %>% 
  inner_join(host_df) %>% 
  select(forest_cover,p_value,patch_no, Land_use, Name, species, z)

#Build Local networks
# Get the names of host and parasite species
Minc <- as.matrix(Minc)

host_species <- rownames(Minc)
parasite_species <- colnames(Minc)

for(p in c(0.01,0.25,0.5)){
  for(f in c(0.1,0.3,0.5,0.7)){

    result_hosts <-  land_h %>% 
      filter(p_value == p) %>% 
      filter(forest_cover == f) 
    
    result_parasites <-land_p %>% 
      filter(p_value == p) %>% 
      filter(forest_cover == f)
    
    # Number of sites
    num_sites <- unique(result_parasites$patch_no)
    
    # List to store the interaction matrices for each site
    var_name <- paste0("ICT_p", p, "_lnet_", f)
    
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


for(p in c(0.01,0.25,0.5)){
  for(f in c(0.1,0.3,0.5,0.7)){
    
    var_name <- paste0("ICT_p", p,"_lnet_", f)
    
    filtered_matrices <- lapply(get(paste0("ICT_p", p,"_lnet_", f)), function(mat) {
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



# Initialize an empty list to store interaction data frames
interaction_list <- list()
idx <- 1  # Index for tracking list position

# Loop through combinations of alpha, p_value, and forest_cover
  for(f in c(0.1, 0.3, 0.5, 0.7)) {
    for(p in c(0.01, 0.25, 0.5)) {
      
      # Get the list of patches
      local_net <- get(paste0("ICT_p", p, "_lnet_", f))
      unique_patch <- names(local_net)
      
      # Loop through each patch matrix
      for(patch in unique_patch) {
        
        # Convert the interaction matrix to a data frame and add metadata
        df <- as.data.frame(as.table(local_net[[patch]])) %>% 
          mutate(patch_id = sub(".*_(.*)", "\\1", patch), 
                 forest_cover = f, 
                 p_value = p) %>% 
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


# Combine all data frames into one
interaction_ict_df <- dplyr::bind_rows(interaction_list)
str(interaction_ict_df)

# Save data frames

write.csv(interaction_ict_df, "./Results_land_new/interaction_ict_df.csv")
