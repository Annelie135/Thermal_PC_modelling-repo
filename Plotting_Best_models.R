d4<-Radial_growth_rate[Radial_growth_rate$Day=="4",]
#choose model
mod = 'pawar_2018'
#get start values
start_vals <- get_start_vals(d4$temperature, d4$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018')
# get limits
low_lims <- get_lower_lims(d4$temperature, d4$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(d4$temperature, d4$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018')

start_vals
low_lims

fit4 <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~pawar_2018(temp = temperature, r_tref, e, eh, topt, tref = 25),
                                    data = d4,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit4
calc_params(fit4) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data4 <- data.frame(temperature = seq(min(d4$temperature), max(d4$temperature), length.out = 100))
preds4<-augment(fit4, newdata = new_data4)

# plot data and model fit
plotpawar<-ggplot(d4, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d4) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds4, col = 'blue') +
  coord_cartesian(ylim = c(0, 0.7))+
  theme_classic(base_size = 14) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm/day)')+
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = "grey", linewidth = 0.6),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.6),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 15), # makes axis titles bigger
        axis.text = element_text(size = 13), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)),   # Increase right margin of y-axis title
        axis.ticks = element_line(colour = "grey"),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())

plot(plotpawar)
##########################
d6<-Radial_growth_rate[Radial_growth_rate$Day=="6",]
mod = 'quadratic_2008'
#get start values
start_vals <- get_start_vals(d6$temperature, d6$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')
# get limits
low_lims <- get_lower_lims(d6$temperature, d6$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')
upper_lims <- get_upper_lims(d6$temperature, d6$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')

start_vals
low_lims

fit6 <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~quadratic_2008(temp = temperature, a, b, c),
                                    data = d6,
                                    iter = c(4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit6
calc_params(fit6) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data6 <- data.frame(temperature = seq(min(d6$temperature), max(d6$temperature), length.out = 100))
preds6<-augment(fit6, newdata = new_data6)

# plot data and model fit
plotquad6_no_title<-ggplot(d6, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d6) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = "grey50") +
  coord_cartesian(ylim = c(0, 0.75))+
  scale_x_continuous(breaks = c(10, 15, 20, 25, 30, 35, 40))+
  theme_classic(base_size = 16) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm/day)')+
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = "grey", linewidth = 0.6),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 16), # makes axis titles bigger
        axis.text = element_text(size = 15), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 15)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 15)),   # Increase right margin of y-axis title
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())



plot(plotquad6_no_title)
########################################
d8<-Radial_growth_rate[Radial_growth_rate$Day=="8",]
#choose model
mod = 'briere2_1999'
#get start values
start_vals <- get_start_vals(d8$temperature, d8$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')
# get limits
low_lims <- get_lower_lims(d8$temperature, d8$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')
upper_lims <- get_upper_lims(d8$temperature, d8$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')

start_vals
low_lims

fit8 <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~briere2_1999(temp = temperature, rmax, ctmax, topt, q10),
                                    data = d8,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit8
calc_params(fit8) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data8 <- data.frame(temperature = seq(min(d8$temperature), max(d8$temperature), length.out = 100))
preds8<-augment(fit8, newdata = new_data8)

# plot data and model fit
Day8brie<-ggplot(d8, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d8) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'grey50') +
  scale_x_continuous(breaks = c(10, 15, 20, 25, 30, 35, 40))+
  coord_cartesian(ylim = c(0, 0.75))+
  theme_classic(base_size = 16) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm/day)')+
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = "grey", linewidth = 0.6),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 16), # makes axis titles bigger
        axis.text = element_text(size = 15), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 15)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 15)),   # Increase right margin of y-axis title
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())



plot(Day8brie)
########################
library(ggplot2)
library(dplyr)
library(broom)

# First, ensure predictions include a label for each model/day
preds_pawar <- augment(fit4, newdata = new_data4) %>% 
  mutate(model = "Pawar 2018", Day = "Day 4")  # example, change if needed

preds_quad6 <- augment(fit6, newdata = new_data6) %>% 
  mutate(model = "Quadratic 2008", Day = "Day 6")

preds_briere8 <- augment(fit8, newdata = new_data8) %>% 
  mutate(model = "Briere2 1999", Day = "Day 8")

# combine predictions
all_preds <- bind_rows(preds_pawar, preds_quad6, preds_briere8)

# combine observed data too (optional)
all_data <- bind_rows(
  d4 %>% mutate(Day = "Day 4"),
  d6 %>% mutate(Day = "Day 6"),
  d8 %>% mutate(Day = "Day 8"))

# plot all together
Growth_preds_all<-ggplot() +
  geom_point(data = all_data, aes(x = temperature, y = Adj_Radial_growth_rate_cm.day, color = Day), size = 2, alpha = 0.7) +
  geom_line(data = all_preds, aes(x = temperature, y = .fitted, color = Day), linewidth = 1) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_x_continuous(breaks = seq(10, 40, 5)) +
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  labs(x = "Temperature (°C)",
       y = "Growth rate (cm/day)")+
  theme_classic(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = c(0.5, 1),   # coordinates inside plot (x = center, y = slightly below top)
        legend.justification = c(0.5, 1),
        legend.direction = "horizontal",
    legend.margin = margin(t = 5), # adds 10 pt top margin
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black"),
    axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
    axis.title.y = element_text(margin = margin(r = 13)))   # Increase right margin of y-axis title))
plot(Growth_preds_all)
ggsave("Growth_preds_allfinal.png", width = 11, height = 7)



############################## DAy 4 weibull ##############################

# plot data and model fit
plotwei<-ggplot(d4, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d4) +
  geom_line(aes(temperature, .fitted), linewidth = 1, d_preds4, col = 'grey50') +
  coord_cartesian(ylim = c(0, 0.75))+
  theme_classic(base_size = 14) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm/day)')+
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = "grey", linewidth = 0.6),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.6),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 15), # makes axis titles bigger
        axis.text = element_text(size = 13), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)),   # Increase right margin of y-axis title
        axis.ticks = element_line(colour = "grey"),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())
plot(plotwei)

#still need to remove background lines

######## Day 6 Briere2 ######
# plot data and model fit
Briere2Day6<-ggplot(d6, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d6) +
  geom_line(aes(temperature, .fitted), linewidth = 1, d_preds6, col = 'grey50') +
  coord_cartesian(ylim = c(0, 0.75))+
  theme_classic(base_size = 14) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm/day)')+
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = "grey", linewidth = 0.6),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.6),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 15), # makes axis titles bigger
        axis.text = element_text(size = 13), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)),   # Increase right margin of y-axis title
        axis.ticks = element_line(colour = "grey"),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())

plot(Briere2Day6)



############### Day 8 quadratic
# plot data and model fit
plot_quad8<-ggplot(d8, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d8) +
  geom_line(aes(temperature, .fitted), linewidth = 1, d_preds8, col = 'grey50') +
  coord_cartesian(ylim = c(0, 0.75))+
  theme_classic(base_size = 14) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm/day)')+
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = "grey", linewidth = 0.6),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.6),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 15), # makes axis titles bigger
        axis.text = element_text(size = 13), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)),   # Increase right margin of y-axis title
        axis.ticks = element_line(colour = "grey"),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())

plot(plot_quad8)


################################################################ FINAL #############
##############################################################################################################
##################################################################################################################
library(ggplot2)
library(dplyr)
library(broom)

# First, ensure predictions include a label for each model/day
preds_weibull <- augment(fit4, newdata = new_data4) %>% 
  mutate(model = "Weibull 1995", Day = "Day 4")  # example, change if needed

preds_briere6 <- augment(fit6, newdata = new_data6) %>% 
  mutate(model = "Briere2 1999", Day = "Day 6")

preds_quad8 <- augment(fit8, newdata = new_data8) %>% 
  mutate(model = "Quadratic 2008", Day = "Day 8")

# combine predictions
all_preds <- bind_rows(preds_weibull, preds_briere6, preds_quad8)

# combine observed data too (optional)
all_data <- bind_rows(
  d4 %>% mutate(Day = "Day 4"),
  d6 %>% mutate(Day = "Day 6"),
  d8 %>% mutate(Day = "Day 8"))

# plot all together
Growth_preds_all<-ggplot() +
  geom_point(data = all_data, aes(x = temperature, y = Adj_Radial_growth_rate_cm.day, color = Day), size = 2, alpha = 0.7) +
  geom_line(data = all_preds, aes(x = temperature, y = .fitted, color = Day), linewidth = 1) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_x_continuous(breaks = seq(10, 40, 5)) +
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  labs(x = "Temperature (°C)",
       y = "Growth rate (cm/day)")+
  theme_classic(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = c(0.5, 1),   # coordinates inside plot (x = center, y = slightly below top)
        legend.justification = c(0.5, 1),
        legend.direction = "horizontal",
        legend.margin = margin(t = 5), # adds 10 pt top margin
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)))   # Increase right margin of y-axis title))
plot(Growth_preds_all)
ggsave("Growth_preds_allfinal25Jan.png", width = 11, height = 7)










