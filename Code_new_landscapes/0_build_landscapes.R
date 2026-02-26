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
                         patch_state=integer(),
                         land= character())

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    l <- read.table(paste0("./Data/landscapes/forestcover_",f,"_pvalue_",p,"_edge.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=1:2500,
                    patch_state=as.vector(l),
                    land= "initial")
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
  dplyr::select(-patch_state) %>% 
  mutate(frag_value = 1 -p_value)

str(land_pv)

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
  )

land_pv<- inner_join(land_pv, final_df)

glimpse(land_pv)

# New Landscape data  ####
landscapes_new <- data.frame(forest_cover=double(),
                         p_value=double(),
                         patch_no=integer(),
                         patch_state=integer(),
                         land= character())

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    l <- read.table(paste0("./Data/landscapes_new/forestcover_",f,"_pvalue_",p,"_edge.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=1:2500,
                    patch_state=as.vector(l),
                    land= "new")
    landscapes_new <- rbind(landscapes_new, l)
  }
}

landscapes_new <- landscapes_new %>%
  dplyr::rename(patch_state="V1") %>%
  left_join(., data.frame(patch_no=1:2500, expand.grid(Y=50:1, X=1:50)))

land_new = landscapes_new %>% 
  #select(- c(X, Y)) %>% 
  mutate(Land_use = 
           case_when(patch_state  == 5 ~ "Core",
                     patch_state  == 4 ~ "Edge",
                     patch_state  == 3 ~ "Edge",
                     patch_state  == 2 ~ "Edge",
                     patch_state  == 1 ~ "Edge",
                     patch_state  ==  0 ~ "Disturbed",
                     TRUE ~ NA)) %>% 
  dplyr::select(-patch_state) %>% 
  mutate(frag_value = 1 -p_value)

str(landscapes_new)
unique(landscapes_new$patch_state)

grouped <- land_new %>%
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
  )

land_new<- inner_join(land_new, final_df)

glimpse(land_new)

land_pv <- full_join(land_pv, land_new)


write.csv(land_pv, "./Data/land_pv.csv")

Land_plot=land_pv %>% 
  filter(land=="new") %>% 
  ggplot(aes(x=X, y=Y, fill=as.factor(Land_use))) +
  geom_tile() +
  facet_grid(p_value~forest_cover) +
  coord_fixed(ratio=1) +
  scale_x_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  scale_y_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  scale_fill_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank())
Land_plot

ggsave("./Results_new/Land_plot_new.png", width = 15, height = 15, units = "cm",  dpi = 600)
