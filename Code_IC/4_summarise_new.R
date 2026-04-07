###################################
# date: 02-2026
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Objective: Analyze different landscape configurations
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
###################################

# Load Packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)

########## Load data ###############
land_pv<- read.csv("./Data/land_pv.csv") %>% column_to_rownames(var = "X.1")

host_info <- read.csv("./Results/host_info.csv") %>% column_to_rownames(var = "X")

prst_info <- read.csv("./Results/prst_info.csv") %>% column_to_rownames(var = "X")

Minc <- read.csv("./Data/M_inc.csv") %>% column_to_rownames(var = "X") %>% 
  rename('juquitiba-like_virus_strain_on576' = 'juquitiba.like_virus_strain_on576')

m_p_traitmatch_ic <- read.csv("./Results_land_new/m_p_traitmatch_alpha.csv") %>% column_to_rownames(var = "X")

land_net_mean <- read.csv("./Results_land_new/alpha_net_mean.csv") %>% column_to_rownames(var = "X")

interaction_land_df <- read.csv("./Results_land_new/interaction_alpha_df.csv") %>% column_to_rownames(var = "X")

# prepare species data 
host_df <- host_info %>%
  rownames_to_column(var= "species") %>% 
  mutate(guild = "resources")

prst_df <- prst_info %>% 
  mutate(guild = "consumers",
         ZoonoticStatus = case_when(
           ZoonoticStatus == 0 ~ "Non-zoonotic",
           ZoonoticStatus == 1 ~ "Zoonotic"))

str(prst_df)

# gather model data
load("./Results_land_new/Sensitivity_ic.RData")


# Calculate parasite global extinction

Ext_new <-  Sensitivity_ic %>%
  filter(guild == "consumers") %>% 
  inner_join(prst_df) %>% 
  select(forest_cover,p_value, Name, species) %>% 
  mutate(Name = as.factor(Name)) %>%  
  group_by(forest_cover, p_value) %>%
  reframe(Prst_richness= as.numeric(n_distinct(Name)),
          Ext_rate = (103 - Prst_richness)/103) %>% 
  mutate(Landscape = paste(forest_cover, p_value, sep = "_")) %>% 
  mutate(frag_value = 1 -p_value) %>% 
  glimpse()
  
# Summarise host data
glimpse(land_pv)

res_h_land <- land_pv %>%
  select(forest_cover, p_value, patch_no, Land_use, core_patch_count, edge_patch_count, disturbed_patch_count) %>% 
  full_join(Sensitivity_ic) %>%
  filter(guild == "resources") %>% 
  mutate(species = as.character(species)) %>% 
  inner_join(host_df) %>% 
  select(forest_cover, p_value, Land_use, core_patch_count, edge_patch_count, disturbed_patch_count, Name, z, pf) %>% 
  dplyr::mutate(Name = as.factor(Name)) %>%  
  group_by(forest_cover, p_value, Land_use, core_patch_count, edge_patch_count, disturbed_patch_count, Name, pf) %>%
  dplyr::summarise(n= n()) %>%
  dplyr::mutate(Abd_type= case_when(
    Land_use == "Core" ~ n/core_patch_count,
    Land_use == "Edge" ~ n/edge_patch_count,
    Land_use == "Disturbed" ~ n/disturbed_patch_count)) %>% 
  ungroup() %>% 
  dplyr::select(- c(core_patch_count, edge_patch_count, disturbed_patch_count)) %>%
  dplyr::rename(host_species = Name) %>% 
  mutate_if(is.character, as.factor) %>% 
  glimpse()

h_abd_land <- land_pv %>% 
  dplyr::rename(patch_id = patch_no) %>% 
  mutate(patch_id = as.character(patch_id)) %>% 
  right_join(interaction_ic_df) %>% 
  unite("Interaction", parasite_species, host_species, sep= "/.") %>% 
  group_by(forest_cover,p_value, Land_use, core_patch_count, edge_patch_count,  disturbed_patch_count, Interaction, value) %>% 
  dplyr::summarise(Count = sum(value), .groups = "drop") %>% 
  separate(Interaction, into = c("parasite_species", "host_species"), sep= "/.") %>% 
  inner_join(res_h_land)  %>% 
  select(forest_cover,p_value, Land_use, parasite_species, Abd_type) %>% 
  group_by(forest_cover,p_value, Land_use, parasite_species) %>% 
  dplyr::summarise(Host_abd = mean(Abd_type), Host_abd_sd = sd(Abd_type)) %>%
  mutate(CoefVar_host = Host_abd_sd/Host_abd) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.factor) %>% 
  glimpse()

#### Parasite and Host attributes 
prst_net <- prst_df %>% 
  select(Name, ZoonoticStatus, Degree, ParGroup) %>% 
  dplyr::rename(parasite_species = Name)

host_net<- host_df %>% 
  dplyr::rename(host_species = Name)

# calculate dilution effect and gather all other outcomes

nh_df_a <- land_pv %>%
  dplyr::rename(patch_id = patch_no) %>% 
  mutate(patch_id = as.character(patch_id)) %>% 
  right_join(interaction_ic_df) %>%
  select(forest_cover,p_value, Land_use, patch_id, parasite_species) %>% 
  group_by(forest_cover,p_value, Land_use, patch_id, parasite_species) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  inner_join(prst_net) %>% 
  mutate(Degree_rel = (Count/Degree)) %>% 
  group_by(forest_cover,p_value, Land_use, parasite_species, ZoonoticStatus, ParGroup) %>%
  dplyr::summarise(Mean_count = mean(Count), SD_rel= sd(Degree_rel), Mean_rel = mean(Degree_rel), .groups = "drop") %>% 
  mutate(zs_dif = Mean_rel/SD_rel) %>% 
  glimpse()

# PCA axis 

pca_data <- new_net_mean_ic %>% 
  select(Mean_size, Mean_connectance, Mean_NODF, Mean_mod) %>%
  scale()  %>% 
  as.data.frame() %>% 
  glimpse()

pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
biplot(pca_result)

axis_pca <- data.frame(pca_result[["rotation"]])

new_net_mean_ic$PCA_axis1 <- pca_result$x[, 1] #connectance positive, size and modularity negative


outcomes_ic <- nh_df_a %>%
  select(- Mean_count) %>%
  inner_join(h_abd_land) %>%
  mutate_if(is.character, as.factor) %>% 
  inner_join(m_p_traitmatch_ic) %>% 
  select(-c(A_mtheta, A_nopartners, Count)) %>% 
  inner_join(Ext_new) %>% 
  inner_join(new_net_mean_ic)


# Analyse outcomes distribution

outcomes_ic$Mean_rel_sc <- as.numeric(scale(outcomes_ic$Mean_rel))
outcomes_ic$Amatch_sc <- as.numeric(scale(outcomes_ic$A_match))
outcomes_ic$Host_abd_sc <- as.numeric(scale(outcomes_ic$Host_abd))
outcomes_ic$Ext_rate_sc <- as.numeric(scale(outcomes_ic$Ext_rate))

hist(outcomes_ic$Mean_rel_sc)
hist(outcomes_ic$Amatch_sc)
hist(outcomes_ic$Host_abd_sc)
hist(outcomes_ic$Ext_rate_sc)

# save data
write.csv(outcomes_ic, "./Results_land_new/outcomes_ic.csv")
