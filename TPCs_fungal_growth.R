library(ggplot2)
library(lme4)
library(performance)
library(car)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(DHARMa)
library(emmeans)
library(broom)
library(MuMIn)
library(nls.multstart)
####################### raw colony growth rate data against temperature  ########################

Radial_growth_ex <- Radial_growth_ex %>%
  filter(!Day %in% c(0))

Radial_growth_ex <- Radial_growth_ex %>%
  filter(!Adj_Radial_growth_rate_cm.day %in% c(NA))


# Colony growth rate across days
Growth_rate<-ggplot(data = Radial_growth_ex, 
                aes(x = temperature, y = Adj_Radial_growth_rate_cm.day, fill = factor(temperature))) +
  geom_bar(stat = "summary", fun = "mean", width = 2, color = "black",linewidth = 0.2 ) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 1, size=unit(0.3, "cm")) +
  labs(x = "Temperature (°C)",y = "Growth rate (cm/day)",fill = "black") +
  facet_wrap(~Day)+
  scale_fill_viridis_d()+
  scale_x_continuous(expand = c(0, 0),breaks = c(10, 15, 20, 25, 30, 35, 40)) +
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 13)+
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = "grey", linewidth = 0.5),
        panel.spacing.x = unit(0.5, "lines"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey95", linewidth =0.2),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.6),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none", legend.title=element_blank()) 
plot(Growth_rate) 


################################################################################################################

library(purrr)
library(rTPC)

#TPC 
d<-Radial_growth_ex[Radial_growth_ex$Day=="4",]
day4_loess<-ggplot(d, aes(x=temperature,y=Adj_Radial_growth_rate_cm.day))+geom_point()+geom_smooth(method = "loess", se=FALSE)
plot(day4_loess)
ggsave("day4_loess.png", width = 13, height = 7)

d_fits <- nest(d, data = c(temperature, Adj_Radial_growth_rate_cm.day)) %>%
  mutate(briere2 = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~briere2_1999(temp = temperature, tmin, tmax, a,b),
                                            data = d,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999') - 5,
                                            start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999') + 5,
                                            lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999'),
                                            upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~gaussian_1987(temp = temperature, rmax, topt, a),
                                             data = d,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         modifiedgaussian = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~modifiedgaussian_2006(temp = temperature, rmax, topt, a, b),
                                                     data = d,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'modifiedgaussian_2006') - 5,
                                                     start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'modifiedgaussian_2006') + 5,
                                                     lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'modifiedgaussian_2006'),
                                                     upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'modifiedgaussian_2006'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         lactin2 = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~lactin2_1995(temp = temperature, a, b, tmax, delta_t),
                                            data = d,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'lactin2_1995') - 10,
                                            start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'lactin2_1995') + 10,
                                            lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'lactin2_1995'),
                                            upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'lactin2_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
          ratkowsky = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~ratkowsky_1983(temp = temperature, tmin, tmax, a, b),
                                              data = d,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'ratkowsky_1983') - 5,
                                              start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'ratkowsky_1983') + 5,
                                              lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'ratkowsky_1983'),
                                              upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~sharpeschoolhigh_1981(temp = temperature, r_tref,e,eh,th, tref = 25),
                                                     data = d,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~quadratic_2008(temp = temperature, a, b, c),
                                              data = d,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008') - 5,
                                              start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008') + 5,
                                              lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~oneill_1972(temp = temperature, rmax, ctmax, topt, q10),
                                           data = d,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'oneill_1972') - 10,
                                           start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'oneill_1972') + 10,
                                           lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'oneill_1972'),
                                           upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'oneill_1972'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         pawar = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~pawar_2018(temp = temperature, r_tref, e, eh, topt, tref = 25),
                                          data = d,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018') - 10,
                                          start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018') + 10,
                                          lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018'),
                                          upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(Adj_Radial_growth_rate_cm.day~weibull_1995(temp = temperature, a,topt,b,c),
                                            data = d,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995'),
                                            supp_errors = 'Y')))

names(d_fits)
glimpse(select(d_fits, 1:7))
d_fits$briere2[[1]]

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', briere2:weibull)

# get predictions using augment
newdata <- tibble(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
d_preds <- d_stack %>%
  mutate(d_stack, preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)
#
d_preds$model_name
# plot
d_predsplotDay414Aug<-ggplot(d_preds, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d) +
  geom_line(aes(temperature, .fitted), col = "blue") +
  facet_wrap(~model_name, labeller = labeller(model_name= label_both), scales = 'free', ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm/day)',
       geom_hline(aes(yintercept = 0), linetype = 2),
       coord_cartesian(ylim = c(0, 0.2)))
plot(d_predsplotDay414Aug)

# plot
g8<-ggplot(d_preds, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d) +
  geom_line(aes(temperature, .fitted), col = "blue") +
  facet_wrap(~model_name, labeller = labeller(model_name= label_both), scales = 'free', ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate (cm²/day)',
       geom_hline(aes(yintercept = 0), linetype = 2),
       coord_cartesian(ylim = c(0, 0.2)))

plot(g8)
ggsave("g8.png", width = 15, height = 7)


AIC_values <- sapply(d_fits, function(model) {
  if (inherits(model, "list") && length(model) > 0) {  
    model <- model[[1]]  # Extract the first sublist (which is the actual `nls` object)
  }
  if (inherits(model, "nls")) {  
    AIC(model)
  } else {
    NA  # Assign NA if it's not a valid model
  }
})
view(AIC_values)
best_model <- names(which.min(AIC_values))
print(best_model)


d_compare <- d_stack %>%
  mutate(
    AIC = map_dbl(fit, ~ if (inherits(.x, "nls")) AIC(.x) else NA_real_),
    BIC = map_dbl(fit, ~ if (inherits(.x, "nls")) BIC(.x) else NA_real_)
  ) %>%
  select(ID, model_name, AIC, BIC)


view(d_compare)



#######################################################  BRIERE 2 #####################################################
d<-Radial_growth_ex[Radial_growth_ex$Day=="4",]
#choose model
mod = 'briere2_1999'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~briere2_1999(temp = temperature, rmax, ctmax, topt, q10),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
Day8brie<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
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

ggsave("Day8briegrey.png", width = 10, height = 6)
# Predictions
d$pred <- predict(fit)
#################
# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)

# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")


model_stats1 <- data.frame(
  model = "briere_1999",
  r_squared = 0.9628,
  adj_r_squared = 0.9415,
  AIC = -31.25,
  AICc = -21.25,
  BIC = -28.83)
################################################  GAUSSIAN   #########################################################

mod = 'gaussian_1987'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'gaussian_1987')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'gaussian_1987')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'gaussian_1987')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~gaussian_1987(temp = temperature, rmax, topt, a),
                                    data = d,
                                    iter = c(4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
Day8gau<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 1))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(Day8gau)

ggsave("Day8gau.png", width = 10, height = 7)


# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)

# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")

model_stats2 <- data.frame(
  model = "gaussian_1987",
  r_squared = 0.892,
  adj_r_squared = 0.8515,
  AIC = -20.47 ,
  AICc = -14.76,
  BIC = -18.53)
####################################################   MODIFIED GAUSSIAN    ############################################
mod = 'modifiedgaussian_2006'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'modifiedgaussian_2006')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'modifiedgaussian_2006')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'modifiedgaussian_2006')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~modifiedgaussian_2006(temp = temperature, rmax, topt, a, b),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
Day8modg<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 1))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(Day8modg)

ggsave("Day6modg.png", width = 10, height = 7)
# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)

# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")

model_stats3 <- data.frame(
  model = "modifiedgaussian_2006",
  r_squared = 0.9488 ,
  adj_r_squared = 0.9195,
  AIC = -27.43 ,
  AICc = -17.43,
  BIC = -25.00)
#####################################################   QUADRATIC   ###################################################
d6<-Radial_growth_rate[Radial_growth_rate$Day=="6",]
mod = 'quadratic_2008'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~quadratic_2008(temp = temperature, a, b, c),
                                    data = d,
                                    iter = c(4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

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
ggsave("plotquad6grey.png", width = 10, height = 6)
# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)

# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")

model_stats4 <- data.frame(
  model = "quadratic_2008",
  r_squared = 0.9518 ,
  adj_r_squared = 0.9337,
  AIC = -30.14 ,
  AICc = -24.43,
  BIC = -28.20)
################################################    SHARPESCHOOLHIGH    #####################################################

#choose model
mod = 'sharpeschoolhigh_1981'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'sharpeschoolhigh_1981')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'sharpeschoolhigh_1981')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~sharpeschoolhigh_1981(temp = temperature, r_tref,e,eh,th, tref = 25),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
plotsharpe<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 1))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(plotsharpe)

# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)
# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")

model_stats5 <- data.frame(
  model = "sharpeschoolhigh_1981",
  r_squared = 0.8719 ,
  adj_r_squared = 0.7987,
  AIC = -16.42,
  AICc = -6.42,
  BIC = -14.00)
####################################################  LACTIN2  ######################################

#choose model
mod = 'lactin2_1995'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'lactin2_1995')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'lactin2_1995')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'lactin2_1995')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~lactin2_1995(temp = temperature, a, b, tmax, delta_t),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
plotlactin<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 1))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(plotlactin)

# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)

# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")

model_stats6 <- data.frame(
  model = "lactin2_1995",
  r_squared = 0.9131,
  adj_r_squared = 0.8635,
  AIC = -21.08 ,
  AICc = -11.08,
  BIC = -18.66)
####################################################    PAWAR   #######################################################

#choose model
mod = 'pawar_2018'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'pawar_2018')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~pawar_2018(temp = temperature, r_tref, e, eh, topt, tref = 25),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
plotpawar<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
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


ggsave("plotpawar.png", width = 10, height = 7)

library(MuMIn)

# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)
# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")
model_stats7 <- data.frame(
  model = "pawar_2018",
  r_squared = 0.8719,
  adj_r_squared = 0.7987,
  AIC = -16.42 ,
  AICc = -6.42,
  BIC = -14.00)

######################################################    RATKOWSKY   ############################################

mod = 'ratkowsky_1983'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'ratkowsky_1983')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'ratkowsky_1983')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'ratkowsky_1983')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~ratkowsky_1983(temp = temperature, tmin, tmax, a, b),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
Day8rat<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 0.75))+
  theme_bw(base_size = 14) +
  labs(x = 'Temperature (ºC)',
       y = 'Inhibition rate (cm²/day)')+
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(Day8rat)

ggsave("Day6rat.png", width = 10, height = 7)
# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)
# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")

model_stats8 <- data.frame(
  model = "ratkowsky_1983",
  r_squared = 0.9482,
  adj_r_squared = 0.9186,
  AIC = -27.28,
  AICc = -17.28,
  BIC = -24.86)

#########################################################   ONEILL   ####################################################

#choose model
mod = 'oneill_1972'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'oneill_1972')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'oneill_1972')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'oneill_1972')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~oneill_1972(temp = temperature, rmax, ctmax, topt, q10),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
Day8one<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 1))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(Day8one)

ggsave("Day8one.png", width = 10, height = 7)

# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)
# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")
model_stats9 <- data.frame(
  model = "oneill_1972",
  r_squared = 0.9088,
  adj_r_squared = 0.8567,
  AIC = -20.5,
  AICc = -10.5,
  BIC = -18.08)

#####################################################   WEIBULL     ############################################

d<-Radial_growth_rate[Radial_growth_rate$Day=="4",]
mod = 'weibull_1995'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~weibull_1995(temp = temperature, a,topt,b,c),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
preds<-augment(fit, newdata = new_data)

# plot data and model fit
plotwei<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'grey50') +
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
        axis.ticks = element_line(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())

plot(plotwei)

ggsave("plotwei_grey.png", width = 10, height = 6)



ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'grey50') +
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
        axis.ticks = element_line(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.spacing.x = unit(1, "lines"),
        legend.title = element_blank())


ggsave("plotpawar_grey.png", width = 10, height = 6)



# Predictions
d$pred <- predict(fit)

# R² and Adjusted R² (Kvalseth, 1985)
rss <- sum((d$Adj_Radial_growth_rate_cm.day - d$pred)^2)
tss <- sum((d$Adj_Radial_growth_rate_cm.day - mean(d$Adj_Radial_growth_rate_cm.day))^2)
r2 <- 1 - (rss / tss)

n <- nrow(d)
p <- length(coef(fit))
adj_r2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - r2)

# AIC and AICc
aic_val <- AIC(fit)
aicc_val <- AICc(fit)
bic_val<-BIC(fit)
# Print metrics
cat("R²:", round(r2, 4), "\n")
cat("Adjusted R² (Kvalseth):", round(adj_r2, 4), "\n")
cat("AIC:", round(aic_val, 2), "\n")
cat("AICc:", round(aicc_val, 2), "\n")
cat("BIC:", round(bic_val, 2), "\n")
model_stats10 <- data.frame(
  model = "weibull_1995",
  r_squared = 0.9058 ,
  adj_r_squared = 0.852,
  AIC = -20.11,
  AICc = -10.11,
  BIC = -17.69)

library(dplyr)
model_statsday8test <- bind_rows(model_stats1, model_stats2, model_stats3, model_stats4, model_stats5, model_stats6, model_stats7, model_stats8, model_stats9,model_stats10)
write.csv(model_statsday8test, "model_statsday8test.csv", row.names = FALSE)
#################################################################################################################################
library(patchwork)

plotwei+plotquad6+Day8brie

#################################################################################################################################
################################################  delta  BIC  Day 4  ############################################################
library(dplyr)
#upload Topt dataset with BIC, AIC and AICc values

# Df= After_4_days

# Calculate delta BIC (ΔBIC)
TPC_Day4 <- TPC_Day4 %>%
  mutate(
    delta_BIC = BIC - min(BIC),
    BIC_weight = exp(-0.5 * delta_BIC) / sum(exp(-0.5 * delta_BIC)))
# Sort from best (lowest BIC)
TPC_Day4 <- TPC_Day4[order(TPC_Day4$BIC), ]
# View
print(TPC_Day4)
################################################################################################################################
################################################  delta  BIC  Day 6  ###########################################################
# Df= TPC_Day6

# Calculate delta BIC (ΔBIC)
TPC_Day6 <- TPC_Day6  %>%
  mutate(
    delta_BIC = BIC - min(BIC),
    BIC_weight = exp(-0.5 * delta_BIC) / sum(exp(-0.5 * delta_BIC)))
# Sort from best (lowest BIC)
TPC_Day6  <- TPC_Day6 [order(TPC_Day6 $BIC), ]
# View
print(TPC_Day6 )

###############################################################################################################################
################################################  delta  BIC  Day 8  ##########################################################
# DF= After_8_days

# Calculate delta BIC (ΔBIC)
TPC_Day8  <- TPC_Day8 %>%
  mutate(
    delta_BIC = BIC - min(BIC),
    BIC_weight = exp(-0.5 * delta_BIC) / sum(exp(-0.5 * delta_BIC)))
# Sort from best (lowest BIC)
TPC_Day8<- TPC_Day8[order(TPC_Day8$BIC), ]
# View
print(TPC_Day8)
##############################################################################################################################
################################################    AIC with correction Day 4  ###############################################
# Calculate delta AICc (ΔAICc)
TPC_Day4 <- TPC_Day4 %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc)))
# Sort from best (lowest AICc)
TPC_Day4 <- TPC_Day4[order(TPC_Day4$AICc), ]
# View
print(TPC_Day4)

##############################################################################################################################
################################################    AIC with correction Day 6  ###############################################
# Df= After_6_days

# Calculate delta AICc (ΔAICc)
TPC_Day6 <- TPC_Day6  %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc)))
# Sort from best (lowest AICc)
TPC_Day6  <- TPC_Day6 [order(TPC_Day6 $AICc), ]
# View
print(TPC_Day6)


#############################################################################################################################
################################################    AIC with correction Day 8  ##############################################

# Dataframe called = After_8_days, renamed to TPC_Day8

# Calculate delta AICc (ΔAICc)
TPC_Day8  <- TPC_Day8 %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc)))
# Sort from best (lowest AICc)
TPC_Day8<- TPC_Day8[order(TPC_Day8$AICc), ]
# View
print(TPC_Day8)

############################################################################################################################
################################################    AIC Day 4  #############################################################

# Calculate delta AIC (ΔAIC)
TPC_Day4 <- TPC_Day4 %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)))
# Sort from best (lowest AIC)
TPC_Day4 <- TPC_Day4[order(TPC_Day4$AIC), ]
# View
print(TPC_Day4)

############################################################################################################################
################################################    AIC Day 6  #############################################################

#Subset Day 6

# Calculate delta AIC (AIC)
TPC_Day6 <- TPC_Day6 %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)))
# Sort from best (lowest AIC)
TPC_Day6  <- TPC_Day6 [order(TPC_Day6 $AIC), ]
# View
print(TPC_Day6)


############################################################################################################################
################################################    AIC Day 8  #############################################################

#Subset Day 8

# Calculate delta AIC (ΔAIC)
TPC_Day8 <- TPC_Day8 %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)))
# Sort from best (lowest BIC)
TPC_Day8 <- TPC_Day8[order(TPC_Day8$AIC), ]
# View
print(TPC_Day8)


############################################################################################################################
library(minpack.lm)
library(nlstools)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)       # for Boot
library(forcats)
library(stringr)
################################################## PLOT BEST MODELS ################################################

d<-Radial_growth_ex[Radial_growth_ex$Day=="4",]
#best model day 4 closest Topt was Sharpeschoolhigh/Pawar but chose Weibull
mod = 'weibull_1995'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~weibull_1995(temp = temperature, a,topt,b,c),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
d_preds<-augment(fit, newdata = new_data)

# plot data and model fit
plotwei<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, d_preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 0.75))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(plotwei)


# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(Adj_Radial_growth_rate_cm.day~weibull_1995(temp = temperature, a,topt,b,c),
                               data = d,
                               start = as.list(coef(fit)),
                               lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995'),
                               upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'weibull_1995'),
                               weights = rep(1, times = nrow(d)))


# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = weibull_1995(temp = temperature, a,topt,b,c))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temperature) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temperature, .fitted), d_preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 0.8))+
  geom_ribbon(aes(temperature, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 13) +
  scale_x_continuous(minor_breaks= c(15,25,35), breaks=c(10,15,20,25,30,35,40), labels = c("10", "15", "20", "25", "30", "35", "40"))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Radial growth rate (cm/day)',
       title = 'Bootstrapped confidence intervals: Radial growth rate \n across temperatures')
p1


# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temperature, .fitted), d_preds, col = 'blue') +
  geom_line(aes(temperature, pred, group = iter), boot1_preds, col = 'blue', alpha = 0.007) +
  coord_cartesian(ylim = c(0, 0.8))+
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 13) +
  scale_x_continuous(minor_breaks= c(15,25,35), breaks=c(10,15,20,25,30,35,40), labels = c("10", "15", "20", "25", "30", "35", "40"))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Radial growth rate (cm/day)',
       title = 'Bootstrapped predictions: Radial growth rate \n across temperatures')

p2
library(patchwork)
p1 + p2

# bootstrap confidence intervals for the extra parameters
extra_params <- calc_params(fit_nlsLM) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate') #extracts the model’s parameter estimates (e.g., r_tref, e, eh, th) from your nlsLM fit.

ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 1000, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params <- left_join(ci_extra_params, extra_params)
#> Joining, by = "param"

ggplot(ci_extra_params, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       subtitle = 'Fusarium TPC; using case resampling')

ci_extra_params
write.csv(ci_extra_params, "CI_day4.csv", row.names = FALSE)

#############################################################################################################################
d<-Radial_growth_ex[Radial_growth_ex$Day=="6",]
#choose model
mod = 'briere2_1999'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999')

start_vals
low_lims
upper_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~briere2_1999(temp = temperature, rmax, ctmax, topt, q10),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
d_preds<-augment(fit, newdata = new_data)

# plot data and model fit
Day6brie<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 1))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(Day6brie)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(Adj_Radial_growth_rate_cm.day~briere2_1999(temp = temperature, rmax, ctmax, topt, q10),
                               data = d,
                               start = as.list(coef(fit)),
                               lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999'),
                               upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'briere2_1999'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = briere2_1999(temp = temperature, rmax, ctmax, topt, q10))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temperature) %>%
  summarise(conf_lower = quantile(pred, 0.025, na.rm = TRUE),
            conf_upper = quantile(pred, 0.975, na.rm = TRUE)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temperature, .fitted), d_preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 0.8))+
  geom_ribbon(aes(temperature, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 13) +
  scale_x_continuous(minor_breaks= c(15,25,35), breaks=c(10,15,20,25,30,35,40), labels = c("10", "15", "20", "25", "30", "35", "40"))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Radial growth rate (cm/day)',
       title = 'Bootstrapped confidence intervals: Radial growth rate \n across temperatures')
p1


# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temperature, .fitted), d_preds, col = 'blue') +
  geom_line(aes(temperature, pred, group = iter), boot1_preds, col = 'blue', alpha = 0.007) +
  coord_cartesian(ylim = c(0, 0.8))+
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 13) +
  scale_x_continuous(minor_breaks= c(15,25,35), breaks=c(10,15,20,25,30,35,40), labels = c("10", "15", "20", "25", "30", "35", "40"))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Radial growth rate (cm/day)',
       title = 'Bootstrapped predictions: Radial growth rate \n across temperatures')

p2
library(patchwork)
p1 + p2

####### 1) Getting the model parameters  

# get parameters of fitted model
param_bact <- broom::tidy(fit_nlsLM) %>%
  select(param = term, estimate)

###### 2) Confidence intervals for model parameters  
# a) Asymptotic CI - Uses the standard error and assumes normal distribution
ci_bact1 <- nlstools::confint2(fit_nlsLM, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')

# b) Profile likelihood  CI - More accurate but slower
ci_bact2 <- nlstools::confint2(fit_nlsLM, trace = TRUE, step = 0.1, maxiter = 1000, tol = 1e-4) %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'profile')
#> Waiting for profiling to be done...

# c) Case bootsrap CI - Resample the rows of your dataset and refit
ci_bact3 <- confint(boot1, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

# d) Residual bootsrap CI - re-sample individuals instead of rows
ci_bact4 <- Boot(fit_nlsLM, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_bact <- bind_rows(ci_bact1, ci_bact2, ci_bact3, ci_bact4) %>%
  left_join(., param_bact)
#> Joining, by = "param"
ci_bact


# plot
ggplot(ci_bact, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_hline(aes(yintercept = conf_lower), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_hline(aes(yintercept = conf_upper), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'Dashed lines are CI of profiling method')


# bootstrap confidence intervals for the extra parameters
extra_params <- calc_params(fit_nlsLM) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate') #extracts the model’s parameter estimates (e.g., r_tref, e, eh, th) from your nlsLM fit.

ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 1000, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params <- left_join(ci_extra_params, extra_params)
#> Joining, by = "param"

ggplot(ci_extra_params, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       subtitle = 'Fusarium TPC; using case resampling')

ci_extra_params
write.csv(ci_extra_params, "CI_day6.csv", row.names = FALSE)



###################################################################################################################################
# Day 8 quadratic 
d<-Radial_growth_ex[Radial_growth_ex$Day=="8",]
mod = 'quadratic_2008'
#get start values
start_vals <- get_start_vals(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')
# get limits
low_lims <- get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')
upper_lims <- get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008')

start_vals
low_lims

fit <- nls.multstart::nls_multstart(Adj_Radial_growth_rate_cm.day~quadratic_2008(temp = temperature, a, b, c),
                                    data = d,
                                    iter = c(4,4,4),
                                    start_lower = start_vals - 5,
                                    start_upper = start_vals + 5,
                                    lower = low_lims,
                                    upper = upper_lims,
                                    supp_errors = 'Y') 
fit
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))
d_preds<-augment(fit, newdata = new_data)

# plot data and model fit
plotquad8<-ggplot(d, aes(temperature, Adj_Radial_growth_rate_cm.day)) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day),d) +
  geom_line(aes(temperature, .fitted), linewidth = 1, preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 1))+
  theme_bw(base_size = 14) +
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
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth =0.2),
        legend.title = element_blank())

plot(plotquad8)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(Adj_Radial_growth_rate_cm.day~quadratic_2008(temp = temperature, a, b, c),
                               data = d,
                               start = as.list(coef(fit)),
                               lower = get_lower_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008'),
                               upper = get_upper_lims(d$temperature, d$Adj_Radial_growth_rate_cm.day, model_name = 'quadratic_2008'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temperature = seq(min(d$temperature), max(d$temperature), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temperature, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temperature) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temperature, .fitted), d_preds, col = 'blue') +
  coord_cartesian(ylim = c(0, 0.8))+
  geom_ribbon(aes(temperature, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 13) +
  scale_x_continuous(minor_breaks= c(15,25,35), breaks=c(10,15,20,25,30,35,40), labels = c("10", "15", "20", "25", "30", "35", "40"))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Radial growth rate (cm/day)',
       title = 'Bootstrapped confidence intervals: Radial growth rate \n across temperatures')
p1


# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temperature, .fitted), d_preds, col = 'blue') +
  geom_line(aes(temperature, pred, group = iter), boot1_preds, col = 'blue', alpha = 0.007) +
  coord_cartesian(ylim = c(0, 0.8))+
  geom_point(aes(temperature, Adj_Radial_growth_rate_cm.day), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 13) +
  scale_x_continuous(minor_breaks= c(15,25,35), breaks=c(10,15,20,25,30,35,40), labels = c("10", "15", "20", "25", "30", "35", "40"))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Radial growth rate (cm/day)',
       title = 'Bootstrapped predictions: Radial growth rate \n across temperatures')

p2
library(patchwork)
p1 + p2

####### 1) Getting the model parameters  

# get parameters of fitted model
param_bact <- broom::tidy(fit_nlsLM) %>%
  select(param = term, estimate)

###### 2) Confidence intervals for model parameters  
# a) Asymptotic CI - Uses the standard error and assumes normal distribution
ci_bact1 <- nlstools::confint2(fit_nlsLM, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')

# b) Profile likelihood  CI - More accurate but slower
ci_bact2 <- nlstools::confint2(fit_nlsLM, trace = TRUE, step = 0.1, maxiter = 100, tol = 1e-4) %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'profile')
#> Waiting for profiling to be done...

# c) Case bootsrap CI - Resample the rows of your dataset and refit
ci_bact3 <- confint(boot1, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

# d) Residual bootsrap CI - re-sample individuals instead of rows
ci_bact4 <- Boot(fit_nlsLM, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_bact <- bind_rows(ci_bact1, ci_bact2, ci_bact3, ci_bact4) %>%
  left_join(., param_bact)
#> Joining, by = "param"

# plot
ggplot(ci_bact, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_hline(aes(yintercept = conf_lower), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_hline(aes(yintercept = conf_upper), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'Dashed lines are CI of profiling method')


# bootstrap confidence intervals for the extra parameters
extra_params <- calc_params(fit_nlsLM) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate') #extracts the model’s parameter estimates (e.g., r_tref, e, eh, th) from your nlsLM fit.

ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 1000, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params <- left_join(ci_extra_params, extra_params)
#> Joining, by = "param"

ggplot(ci_extra_params, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       subtitle = 'Fusarium TPC; using case resampling')

ci_extra_params
write.csv(ci_extra_params, "CI_day8.csv", row.names = FALSE)

#############################################################################################################################
#############################################################################################################################
#combined CI_day4.csv CI_day6.csv CI_day8.csv into one dataframe called = Params_CI_all 
#import the datasheet Params_CI_all from TPC_Models(stats).xlsx 

# filter only Topt values
topt_data <- Params_CI_all %>% filter(param == "topt")

# plot
Topt_CI<-ggplot(topt_data, aes(x = factor(day), y = estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = conf_lower, ymax = conf_upper), width = 0.2, color = "blue") +
  labs(
    title = "Growth rate Topt Estimates",
    x = "Day",
    y = "Temperature (°C)"
  ) +
  coord_flip()+
  theme_minimal(base_size = 14)+
  theme(plot.background = element_rect(fill = "white", color = NA),
        axis.title = element_text(size = 15), # makes axis titles bigger
        axis.text = element_text(size = 13), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)))   # Increase right margin of y-axis title)
plot(Topt_CI)
ggsave("Topt_CI.png", width = 15, height = 7)

#plot Rmax
# filter only Rmax values
rmax_data <- Params_CI_all %>% filter(param == "rmax")

# plot
Rmax_CI<-ggplot(rmax_data, aes(x = factor(day), y = estimate)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = conf_lower, ymax = conf_upper), width = 0.2, color = "red") +
  labs(
    title = "Growth rate Rmax Estimates",
    x = "Day",
    y = "Rmax: growth rate (cm/day)"
  ) +
  coord_flip()+
  theme_minimal(base_size = 14)+
  theme(plot.background = element_rect(fill = "white", color = NA),
        axis.title = element_text(size = 15), # makes axis titles bigger
        axis.text = element_text(size = 13), # makes tick labels bigger
        axis.title.x = element_text(margin = margin(t = 13)),  # Increase top margin of x-axis title
        axis.title.y = element_text(margin = margin(r = 13)))   # Increase right margin of y-axis title)
plot(Rmax_CI)
ggsave("Rmax_CI.png", width = 15, height = 7)



