####################################
# title: Equilibrium state from ecocoevolutionary model without rewiring and different strength of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# date: 03-2025
###################################

#Load packages

pacman::p_load(tidyverse, ggplot2, ggpubr, rstatix, cowplot, ggeffects, ggpredict, performance)


# Load data

time_data <- read.csv("Output_v2b/dt_metanetwork_fc0.3_p0.5_er1_ec0.25_cr1_cc0.4_mr0.5_mc0.7_phi0.5_alpha0.2_eps5.csv")

str(time_data)

# Calculate abundance
unique(time_data$guild)
time_r <- time_data %>%
  filter(abundance != is.na(abundance)) %>% 
  filter(guild == "resources") %>% 
  mutate(species = as.factor(species)) %>% 
  select(t, species, abundance, z_mean, z_sd, guild) %>% 
  ggplot(aes(x= t, y = abundance, fill= species, color = species))+
  geom_line(alpha = 0.5)+
  labs(title = "Host species",
       x = "Time",
       y = "Mean Abundance") +
  theme_bw() +
  theme(title = element_text(size= 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y =element_blank()
  )

time_r


time_c <- time_data %>%
  filter(abundance != is.na(abundance)) %>% 
  filter(guild == "consumers") %>% 
  mutate(species = as.factor(species)) %>% 
  select(t, species, abundance, z_mean, z_sd, guild) %>% 
  ggplot(aes(x= t, y = abundance, fill= species, color = species))+
  geom_line(alpha = 0.5)+
  labs(title = "Parasite species",
                             x = "Time",
                             y = "Mean Abundance") +
  theme_bw() +
  theme(title = element_text(size= 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y =element_blank()
  )

## Merge plots ####

grid <- ggarrange(time_r, time_c, nrow = 2, ncol = 1, widths = c(1, 3))

grid

png(filename="./Figures/time_grid.png", width=20, height=25, units="cm", res=600)
plot(grid)
dev.off()

# Calculate trait
unique(time_data$guild)

time_r <- time_data %>%
  filter(z_mean != is.na(z_mean)) %>% 
  filter(guild == "resources") %>% 
  mutate(species = as.factor(species)) %>% 
  select(t, species, abundance, z_mean, z_sd, guild) %>% 
  ggplot(aes(x= t, y = z_mean, fill= species, color = species))+
  geom_line(alpha = 0.5)+
  labs(title = "Host species",
       x = "Time",
       y = "Mean trait") +
  theme_bw() +
  theme(title = element_text(size= 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y =element_blank()
  )

time_r


time_c <- time_data %>%
  filter(abundance != is.na(z_mean)) %>% 
  filter(guild == "consumers") %>% 
  mutate(species = as.factor(species)) %>% 
  select(t, species, abundance, z_mean, z_sd, guild) %>% 
  ggplot(aes(x= t, y = z_mean, fill= species, color = species))+
  geom_line(alpha = 0.5)+
  labs(title = "Parasite species",
       x = "Time",
       y = "Mean trait") +
  theme_bw() +
  theme(title = element_text(size= 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y =element_blank()
  )

## Merge plots ####


grid <- ggarrange(time_r, time_c, nrow = 2, ncol = 1, widths = c(1, 3))

grid

png(filename="./Figures/time_trait_grid.png", width=20, height=25, units="cm", res=600)
plot(grid)
dev.off()




time_c
