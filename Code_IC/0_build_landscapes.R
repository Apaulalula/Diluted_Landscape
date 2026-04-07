###################################
# date: 02-2026
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Objective: Analyze different landscape configurations
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
###################################

# Build landscapes

# Initial Landscape data  ####
landscapes <- data.frame(forest_cover=double(),
                         p_value=double(),
                         patch_no=integer(),
                         patch_state=integer()
                         )

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    l <- read.table(paste0("./Data/landscapes_new/forestcover_",f,"_pvalue_",p,"_edge.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=1:2500,
                    patch_state=as.vector(l)
                    )
    landscapes <- rbind(landscapes, l)
  }
}

landscapes <- landscapes %>%
  dplyr::rename(patch_state="V1") %>%
  left_join(., data.frame(patch_no=1:2500, expand.grid(Y=50:1, X=1:50)))

land_pv = landscapes %>% 
  #select(- c(X, Y)) %>% 
  mutate(Land_use = 
           case_when(patch_state  == 5 ~ "Core",
                     patch_state  == 4 ~ "Edge",
                     patch_state  == 3 ~ "Edge",
                     patch_state  == 2 ~ "Edge",
                     patch_state  == 1 ~ "Edge",
                     patch_state  ==  0 ~ "Disturbed",
                     TRUE ~ NA)) %>% 
  dplyr::select(-patch_state) 

grouped <- land_pv %>%
  select(p_value, forest_cover, Land_use, patch_no) %>% 
  group_by(p_value, forest_cover, Land_use) %>%
  dplyr::summarise(count = n(), .groups = 'drop')

# Creating new columns for forest and non-forest patch counts
grouped <- grouped %>%
  mutate(
    core_patch_count = ifelse(Land_use == 'Core', count, 0),
    edge_patch_count = ifelse(Land_use == 'Edge', count, 0),
    disturbed_patch_count = ifelse(Land_use == 'Disturbed', count, 0)
  )

# Summing up the patch counts for each forest_cover and fragmentation_level
final_df <- grouped %>%
  group_by(p_value, forest_cover) %>%
  dplyr::summarise(
    core_patch_count = sum(core_patch_count),
    edge_patch_count = sum(edge_patch_count),
    disturbed_patch_count = sum(disturbed_patch_count),
    .groups = 'drop'
  ) %>% 
  as.data.frame() %>% 
  glimpse()

land_pv<- inner_join(land_pv, final_df) %>% 
  select(-c(X,Y)) %>% 
  glimpse()


