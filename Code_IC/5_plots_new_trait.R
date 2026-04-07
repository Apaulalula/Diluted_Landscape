###################################
# date: 02-2026
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Objective: Analyze different landscape configurations
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
###################################

######## Load Packages ######## 
# Load Packages
pacman::p_load(tidyverse, ggplot2, ggpubr)

######## Load Data ######## 
outcomes_ict <- read.csv("./Results_land_new/outcomes_ict.csv") %>% column_to_rownames(var = "X")

# Work on data set ################################### 

# Verify unique levels of categorical variables
outcomes_ict$parasite_species <- factor(outcomes_ict$parasite_species)
outcomes_ict$Landscape <- factor(outcomes_ict$Landscape)


######## Build Plots ######## 

# Extinction plot ####
str(Ex_plot_data)

Ex_plot_data <- outcomes_ict %>%  
  group_by(forest_cover, p_value) %>%
  dplyr::summarise(
    mean_value = mean(Ext_rate),
    se_value = sd(Ext_rate) / sqrt(n()),
    .groups = "drop"
  )

Ex_plot <- ggplot(Ex_plot_data, aes(x= as.factor(forest_cover), y= mean_value)) + 
  geom_pointrange(aes(ymin = mean_value - se_value, ymax = mean_value + se_value, color= as.factor(p_value))) +
  geom_line(aes(group = as.factor(p_value), color= as.factor(p_value)), linetype= "dashed") +
  labs(title = "Traits from uniform distribution")+
  scale_color_manual(values=c("0.01" = '#C2A878', "0.25" = '#6E633D', "0.5"= '#355834'), name=NULL) +
  labs( x = "Amount of forest cover", y= "Parasite extinction rate")+
  theme_classic()+
  theme(legend.position= "none",#c(0.8,0.9),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
Ex_plot


png(filename="./Results_new/ext_md_ict.png", width=10, height=8, units="cm", res=600)
plot(Ex_plot)
dev.off()

# Network structure plots ####
df_plot <- outcomes_ict %>%
  mutate(
    forest_cover = as.factor(forest_cover),
    Land_use = as.factor(Land_use)
  ) %>%
  group_by(forest_cover, p_value) %>%
  dplyr::summarise(
    mean_value = mean(PCA_axis1),
    se_value = sd(PCA_axis1) / sqrt(n()), .groups = "drop"
  )


PC1_plot <- df_plot %>% 
  ggplot(aes(x= as.factor(forest_cover), y= mean_value)) +
  geom_point(aes(color= as.factor(p_value)),size = 1, alpha= 0.5) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value, color= as.factor(p_value)), width = 0.2, alpha= 0.5, size = 1) +
  geom_line(aes(group = as.factor(p_value), color= as.factor(p_value)),linetype ="dotted") +  # Adding lines for each Land_use group
  labs(title = "Traits from uniform distribution")+
  scale_color_manual(values=c("0.01" = '#C2A878', "0.25" = '#6E633D', "0.5"= '#355834'), name=NULL) +
  labs( x = "Amount of forest cover", y= "Network structure PC1")+
  theme_classic()+
  theme(legend.position= "none",#c(0.8,0.9),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
PC1_plot

df_plot <- outcomes_ict %>%
  mutate(
    forest_cover = as.factor(forest_cover),
    Land_use = as.factor(Land_use)
  ) %>%
  group_by(forest_cover, Land_use) %>%
  dplyr::summarise(
    mean_value = mean(PCA_axis1),
    se_value = sd(PCA_axis1) / sqrt(n()), .groups = "drop"
  ) %>% 
  glimpse()

PC1_plot1 <- ggplot(df_plot, aes(x = forest_cover, y = mean_value)) + 
  geom_point(aes(color = Land_use), size = 1.2, alpha = 0.7, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value, color = Land_use),
                width = 0.2, alpha = 0.6, size = 0.8, position = position_dodge(0.3)) +
  geom_line(aes(group = Land_use, color = Land_use), linetype = "dotted", position = position_dodge(0.3)) +
  labs(title = "Traits from uniform distribution")+
  scale_color_manual(values = c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name = "Patch land use") +
  labs(x = "Amount of forest cover", y = "Network structure PC1") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10)
  )

PC1_plot1

# Combine the grid and the y-axis label
final_plot<- plot_grid(PC1_plot, PC1_plot1, ncol = 1)

# Display the final plot
print(final_plot)

png(filename="./Results_new/PCA_ns_ict.png", width=12, height=20, units="cm", res=600)
plot(final_plot)
dev.off()

# Outcomes Plots #######################################
## Relative degree plots #####
### Plot relation between relative degree,zoonotic status and extinction rate  ####
dg_plot1 <-outcomes_ict %>% 
  select(-p_value) %>% 
  ggplot(aes(x= Ext_rate, y= Mean_rel_sc, color= ZoonoticStatus)) + 
  geom_point(alpha= 0.1)+
  geom_smooth(method = lm, se = TRUE)+
  scale_color_manual(values=c("Zoonotic" = '#FF6F59', "Non-zoonotic" = '#254441'), name=NULL) +
  #facet_wrap(~IC)+
  labs(title = "Traits from uniform distribution")+
  labs( x = "Parasite extinction rate", y= "Relative degree (scaled)")+
  theme_classic()+
  theme(legend.position="none", #c(0.8,0.15),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal", 
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
dg_plot1

png(filename="./Results_new/DG_plot_ict.png", width=15, height=7, units="cm", res=600)
plot(dg_plot1)
dev.off()

### Plot relation between relative degree, land use and host abundance  ####
dg_plot <-outcomes_ict %>% 
  select(-p_value) %>% 
  ggplot(aes(x= Host_abd, y= Mean_rel)) + 
  geom_point(alpha= 0.5, aes(color= as.factor(Land_use)))+
  geom_smooth(method = lm, se = TRUE, color= "grey50")+
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Patch land use")+
  labs(title = "Traits from uniform distribution")+
  # facet_wrap(~IC)+
  labs( x = "Host abundance", y= "Relative degree (scaled)")+
  theme_classic()+
  theme(legend.position="none",
        #legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
dg_plot

# Combine the grid and the y-axis label
final_plot<- plot_grid(dg_plot1, dg_plot, ncol = 1)

# Display the final plot
print(final_plot)

png(filename="./Results_new/P_rd_ict.png", width=12, height=20, units="cm", res=600)
plot(final_plot)
dev.off()

## Trait matching plots ####
### Plot relation between trait matching, forest cover and network structure ####
TM_plot1 <-outcomes_ict %>% 
  select(-p_value) %>% 
  ggplot(aes(x= PCA_axis1, y= Amatch_sc)) + 
  geom_point(alpha= 0.15, aes(color= as.factor(Land_use)))+
  geom_smooth(method = lm, se = TRUE, color= "grey50")+
  #scale_color_continuous_sequential(palette = "Greens", h2 = 40)+
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Patch land use")+
#  scale_color_manual(values=c("0.1" = '#CFE0BC', "0.3" = '#7FA653',"0.5" = '#63783D',"0.7" = '#234D21'), name=NULL)+
  # facet_wrap(~IC)+
  labs(title = "Traits from uniform distribution")+
  labs( x = "Network structure PC1", y= "Trait matching (scaled)")+
  theme_classic()+
  theme(legend.position="none",
        #legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
TM_plot1

### Plot relation between trait matching, habitat type and host abundance ####
TM_plot2 <- outcomes_ict %>% 
  select(-p_value) %>% 
  ggplot(aes(x= Host_abd, y= Amatch_sc)) + 
  geom_point(alpha= 0.15, aes(color= as.factor(Land_use)))+
  geom_smooth(method = "lm", se = TRUE, aes(color= Land_use))+
  #scale_color_continuous_sequential(palette = "Greens", h2 = 40)+
  # facet_wrap(~IC)+
  labs(title = "Traits from uniform distribution")+
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Patch land use")+
  labs( x = "Host relative abundance", y= "Trait matching (scaled)")+
  theme_classic()+
  theme(legend.position="none",
        #legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
TM_plot2

# Combine the grid and the y-axis label
final_plot<- plot_grid(TM_plot1, TM_plot2, ncol = 1)

# Display the final plot
print(final_plot)

png(filename="./Results_new/P_tm_ict.png", width=12, height=20, units="cm", res=600)
plot(final_plot)
dev.off()
