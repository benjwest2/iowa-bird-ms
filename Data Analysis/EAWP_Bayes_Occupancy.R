# Ben West's occupancy model

setwd('C:/Users/benjwest/Box Sync/SIFM Bird Data All Years/bwest-bird-density')

library(rjags)
library(tidyverse)
library(plotrix)

logit <- function(x){ln(x/(1-x))}

#This is for EAWP####
occ_dat_org <- read.csv("Data Wrangling SIFM/EAWP_Occ_2019_Wide.csv", header = T) %>%
  rowwise() %>%
  mutate(tssr_sq_vis_1 = tssr_vis_1^2) %>%
  mutate(tssr_sq_vis_2 = tssr_vis_2^2) %>%
  mutate(jday_sq_vis_1 = jday_vis_1^2) %>%
  mutate(jday_sq_vis_2 = jday_vis_2^2) %>%
  as.data.frame()

# Center and scale variables
variables_scaled <- as.data.frame(scale(occ_dat_org[,c(4, 6:8, 10:22)]))
occ_dat <- cbind(occ_dat_org[,c(1:3, 5, 9)], variables_scaled)

summary(variables_scaled$jday_sq_vis_2)

# Initialize the BUGS model####
occ_model <- jags.model(
  'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_Bayes_Occ.bug',
  data = list(n = nrow(occ_dat),
              num_tree_spp = occ_dat$num_tree_spp,
              live_basal = occ_dat$live_basal,
              oak_prop = occ_dat$oak_prop,
              can_close_prop = occ_dat$can_clos_prop,
              shrub_dens = occ_dat$shrub_dens,
              prop_forest_1km = occ_dat$prop_forest_1km,
              prop_forest_10km = occ_dat$prop_forest_10km,
              tssr_vis_1 = occ_dat$tssr_vis_1,
              jday_vis_1 = occ_dat$jday_vis_1,
              obs_vis_1 = occ_dat$obs_vis_1,
              tssr_vis_2 = occ_dat$tssr_vis_2,
              jday_vis_2 = occ_dat$jday_vis_2,
              obs_vis_2 = occ_dat$obs_vis_2,
              detect_vis_1 = occ_dat$detect_vis_1,
              detect_vis_2 = occ_dat$detect_vis_2,
              tssr_sq_vis_1 = occ_dat$tssr_sq_vis_1,
              tssr_sq_vis_2 = occ_dat$tssr_sq_vis_2,
              jday_sq_vis_1 = occ_dat$jday_sq_vis_1,
              jday_sq_vis_2 = occ_dat$jday_sq_vis_2),
    n.chains=3,
    n.adapt=100)

#Making a data frame that is the key for my different parameters####
#This next vector can also be used in my model fit
b <- c('b00', 'b01','b02', 'b03','b04', 'b05','b06', 'b07', 'b08','b09', 'b10','b11','b12')
param_name <- c("intercept", "tree_spp_richness", "live_basal_area", "oak_prop", "canopy_closure", "shrub_density", "prop_forest_in_1km", "prop_forest_in_10km", "minutes_past_sunrise", "mins_past_sr_squared", 'julian_day','julian_day_squared','observer')
det_or_occ <- c('intercept',rep('occupancy',7),rep('detection',5))

param_key <- cbind(b,param_name,det_or_occ) %>%
  as_tibble()

# sample the MCMC chains 10000 times as burn-in, ignore the samples
update(occ_model, 10000)

# sample another 10000 times, keep the variables of interest####
EAWP.fit <- coda.samples(occ_model, 
                           c(b, "true_occ"),
                           10000,
                           thin=1)

#Get summaries of model####

#Mean occupancy probabilities
occ_probs <- do.call(rbind.data.frame, EAWP.fit)[14:514] %>%
  summarize_all(list(mean)) %>%
  matrix(nrow=1) %>%
  t() %>%
  as_tibble() %>%
  rename(mean_occ_prob = V1) %>%
  bind_cols(occ_dat_org) %>%
  select(point_id, mean_occ_prob) %>%
  mutate(mean_occ_prob = unlist(mean_occ_prob))

#Histogram of mean occupancy probs
ggplot(data = occ_probs, aes(x = mean_occ_prob))+
  geom_histogram(binwidth = 0.1,color="black", fill="gray80")+
  theme_classic()+
  xlab("Mean Occupancy Probability")+
  ylab("Count")+
  scale_x_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))

occ_order <- occ_probs %>%
  arrange(desc(mean_occ_prob))

occ_order[381,]
occ_order[382,]

#Getting means and standard errors
EAWP_stats <- summary(EAWP.fit)$statistics %>%
  as.data.frame() %>%
  rownames_to_column(var = "b")%>%
  as_tibble()%>%
  filter(str_detect(b, "^b")) %>%
  select(-`Time-series SE`, -SD)%>%
  rename(SE_Naive = `Naive SE`) %>%
  rowwise() %>%
  mutate(Mean = round(Mean, 3))%>%
  mutate(SE_Naive = round(SE_Naive,3)) %>%
  left_join(param_key) %>%
  select(param_name, Mean, SE_Naive, det_or_occ, b) %>%
  mutate(abs_mean = abs(Mean)) %>%
  arrange(desc(abs_mean)) %>%
  select(-abs_mean)

#Getting credible intervals, determining their width, and seeing if they cross zero
EAWP_quant <- summary(EAWP.fit)$quantiles %>%
  as.data.frame() %>%
  rownames_to_column(var = "b")%>%
  as_tibble() %>%
  filter(str_detect(b, "^b")) %>%
  select(b, `2.5%`, `97.5%`) %>%
  rowwise() %>%
  mutate(ci_width = `97.5%` - `2.5%`) %>%
  mutate(low_pos =ifelse(`2.5%` > 0, 1, 0))%>%
  mutate(up_pos = ifelse(`97.5%` > 0, 1, 0))%>%
  mutate(ci_cross_zero = ifelse(low_pos == up_pos,1,0))%>%
  arrange(desc(ci_cross_zero), ci_width)%>%
  mutate(ci_cross_zero = as.character(ci_cross_zero))%>%
  mutate(ci_cross_zero = ifelse(ci_cross_zero == "0", "yes","no"))%>%
  left_join(param_key) %>%
  select(param_name, `2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ,b) %>%
  mutate(`2.5%` = round(`2.5%`,3))%>%
  mutate(`97.5%` = round(`97.5%`,3))%>%
  mutate(ci_width = round(ci_width,3))
  
#Putting together means, credible intervals, and associated stats
EAWP_summary <- EAWP_stats %>%
  left_join(EAWP_quant) %>%
  select(param_name, Mean, SE_Naive,`2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ,b)

EAWP_summary

write.csv(occ_probs, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_occ_probs.csv')
write.csv(EAWP_summary, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_model_summary.csv')

#Model diagnostics
gelman.diag(EAWP.fit)


#Going to try KEWA####

occ_dat_org2 <- read.csv("Data Wrangling SIFM/KEWA_Occ_2019_Wide.csv", header = T) %>%
  rowwise() %>%
  mutate(tssr_sq_vis_1 = tssr_vis_1^2) %>%
  mutate(tssr_sq_vis_2 = tssr_vis_2^2) %>%
  mutate(jday_sq_vis_1 = jday_vis_1^2) %>%
  mutate(jday_sq_vis_2 = jday_vis_2^2) %>%
  as.data.frame()

# Center and scale variables
variables_scaled2 <- as.data.frame(scale(occ_dat_org2[,c(4, 6:8, 10:22)]))
occ_dat2 <- cbind(occ_dat_org2[,c(1:3, 5, 9)], variables_scaled2)

head(variables_scaled2)

# Initialize the BUGS model####
occ_model2 <- jags.model(
  'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_Bayes_Occ.bug',
  data = list(n = nrow(occ_dat2),
              num_tree_spp = occ_dat2$num_tree_spp,
              live_basal = occ_dat2$live_basal,
              oak_prop = occ_dat2$oak_prop,
              can_close_prop = occ_dat2$can_clos_prop,
              shrub_dens = occ_dat2$shrub_dens,
              prop_forest_1km = occ_dat2$prop_forest_1km,
              prop_forest_10km = occ_dat2$prop_forest_10km,
              tssr_vis_1 = occ_dat2$tssr_vis_1,
              jday_vis_1 = occ_dat2$jday_vis_1,
              obs_vis_1 = occ_dat2$obs_vis_1,
              tssr_vis_2 = occ_dat2$tssr_vis_2,
              jday_vis_2 = occ_dat2$jday_vis_2,
              obs_vis_2 = occ_dat2$obs_vis_2,
              detect_vis_1 = occ_dat2$detect_vis_1,
              detect_vis_2 = occ_dat2$detect_vis_2,
              tssr_sq_vis_1 = occ_dat2$tssr_sq_vis_1,
              tssr_sq_vis_2 = occ_dat2$tssr_sq_vis_2,
              jday_sq_vis_1 = occ_dat2$jday_sq_vis_1,
              jday_sq_vis_2 = occ_dat2$jday_sq_vis_2),
  n.chains=3,
  n.adapt=100)

update(occ_model2, 10000)
KEWA.fit <- coda.samples(occ_model2, 
                         c(b, "true_occ"),
                         10000,
                         thin=1)


#Get summaries of model####

#Mean occupancy probabilities
occ_probs2 <- do.call(rbind.data.frame, KEWA.fit)[14:514] %>%
  summarize_all(list(mean)) %>%
  matrix(nrow=1) %>%
  t() %>%
  as_tibble() %>%
  rename(mean_occ_prob = V1) %>%
  bind_cols(occ_dat_org2) %>%
  select(point_id, mean_occ_prob) %>%
  mutate(mean_occ_prob = unlist(mean_occ_prob))

#Histogram of mean occupancy probs
ggplot(data = occ_probs2, aes(x = mean_occ_prob))+
  geom_histogram(binwidth = 0.1)+
  theme_classic()

#Getting means and standard errors
KEWA_stats <- summary(KEWA.fit)$statistics %>%
  as.data.frame() %>%
  rownames_to_column(var = "b")%>%
  as_tibble()%>%
  filter(str_detect(b, "^b")) %>%
  select(-`Time-series SE`, -SD)%>%
  rename(SE_Naive = `Naive SE`) %>%
  rowwise() %>%
  mutate(Mean = round(Mean, 3))%>%
  mutate(SE_Naive = round(SE_Naive,3)) %>%
  left_join(param_key) %>%
  select(param_name, Mean, SE_Naive, det_or_occ, b) %>%
  mutate(abs_mean = abs(Mean)) %>%
  arrange(desc(abs_mean)) %>%
  select(-abs_mean)

#Getting credible intervals, determining their width, and seeing if they cross zero
KEWA_quant <- summary(KEWA.fit)$quantiles %>%
  as.data.frame() %>%
  rownames_to_column(var = "b")%>%
  as_tibble() %>%
  filter(str_detect(b, "^b")) %>%
  select(b, `2.5%`, `97.5%`) %>%
  rowwise() %>%
  mutate(ci_width = `97.5%` - `2.5%`) %>%
  mutate(low_pos =ifelse(`2.5%` > 0, 1, 0))%>%
  mutate(up_pos = ifelse(`97.5%` > 0, 1, 0))%>%
  mutate(ci_cross_zero = ifelse(low_pos == up_pos,1,0))%>%
  arrange(desc(ci_cross_zero), ci_width)%>%
  mutate(ci_cross_zero = as.character(ci_cross_zero))%>%
  mutate(ci_cross_zero = ifelse(ci_cross_zero == "0", "yes","no"))%>%
  left_join(param_key) %>%
  select(param_name, `2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ,b) %>%
  mutate(`2.5%` = round(`2.5%`,3))%>%
  mutate(`97.5%` = round(`97.5%`,3))%>%
  mutate(ci_width = round(ci_width,3))

#Putting together means, credible intervals, and associated stats
KEWA_summary <- KEWA_stats %>%
  left_join(KEWA_quant) %>%
  select(param_name, Mean, SE_Naive,`2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ,b)

KEWA_summary

write.csv(occ_probs2, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/KEWA_occ_probs.csv')
write.csv(KEWA_summary, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/KEWA_model_summary.csv')

#KEWA Plot####
occ_probs2 <- read_csv('Data Exploration and Analysis/EAWP Bayesian Occupancy/KEWA_occ_probs.csv')
ggplot(data = occ_probs2, aes(x = mean_occ_prob))+
  geom_histogram(binwidth = 0.1,color="black", fill="gray80")+
  theme_classic()+
  xlab("Mean Occupancy Probability")+
  ylab("Count")+
  scale_x_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))


#EAWP Model; fewer parameters####

b2 <- c('b00', 'b01','b02', 'b03','b04', 'b05','b06')
param_name2 <- c("intercept", "live_basal_area", "oak_prop", "canopy_closure", "prop_forest_in_1km", "prop_forest_in_10km",'observer')
det_or_occ2 <- c('intercept',rep('occupancy',5),'detection')

param_key2 <- cbind(b2,param_name2,det_or_occ2) %>%
  as_tibble()

# Initialize the BUGS model####
occ_model3 <- jags.model(
  'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_Bayes_Occ2.txt',
  data = list(n = nrow(occ_dat),
              live_basal = occ_dat$live_basal,
              oak_prop = occ_dat$oak_prop,
              can_close_prop = occ_dat$can_clos_prop,
              prop_forest_1km = occ_dat$prop_forest_1km,
              prop_forest_10km = occ_dat$prop_forest_10km,
              obs_vis_1 = occ_dat$obs_vis_1,
              obs_vis_2 = occ_dat$obs_vis_2,
              detect_vis_1 = occ_dat$detect_vis_1,
              detect_vis_2 = occ_dat$detect_vis_2),
  n.chains=3,
  n.adapt=100)

update(occ_model3, 10000)
EAWP2.fit <- coda.samples(occ_model3, 
                         c(b2, "true_occ"),
                         10000,
                         thin=1)


#Get summaries of model####

#Mean occupancy probabilities
occ_probs3 <- do.call(rbind.data.frame, EAWP2.fit)[8:508] %>%
  summarize_all(list(mean)) %>%
  matrix(nrow=1) %>%
  t() %>%
  as_tibble() %>%
  rename(mean_occ_prob = V1) %>%
  bind_cols(occ_dat_org2) %>%
  select(point_id, mean_occ_prob) %>%
  mutate(mean_occ_prob = unlist(mean_occ_prob))

#Histogram of mean occupancy probs
ggplot(data = occ_probs3, aes(x = mean_occ_prob))+
  geom_histogram(binwidth = 0.1)+
  theme_classic()

#Getting means and standard errors
EAWP2_stats <- summary(EAWP2.fit)$statistics %>%
  as.data.frame() %>%
  rownames_to_column(var = "b2")%>%
  as_tibble()%>%
  filter(str_detect(b2, "^b")) %>%
  select(-`Time-series SE`, -SD)%>%
  rename(SE_Naive = `Naive SE`) %>%
  rowwise() %>%
  mutate(Mean = round(Mean, 3))%>%
  mutate(SE_Naive = round(SE_Naive,3)) %>%
  left_join(param_key2) %>%
  select(param_name2, Mean, SE_Naive, det_or_occ2, b2) %>%
  mutate(abs_mean = abs(Mean)) %>%
  arrange(desc(abs_mean)) %>%
  select(-abs_mean)

#Getting credible intervals, determining their width, and seeing if they cross zero
EAWP2_quant <- summary(EAWP2.fit)$quantiles %>%
  as.data.frame() %>%
  rownames_to_column(var = "b2")%>%
  as_tibble() %>%
  filter(str_detect(b2, "^b")) %>%
  select(b2, `2.5%`, `97.5%`) %>%
  rowwise() %>%
  mutate(ci_width = `97.5%` - `2.5%`) %>%
  mutate(low_pos =ifelse(`2.5%` > 0, 1, 0))%>%
  mutate(up_pos = ifelse(`97.5%` > 0, 1, 0))%>%
  mutate(ci_cross_zero = ifelse(low_pos == up_pos,1,0))%>%
  arrange(desc(ci_cross_zero), ci_width)%>%
  mutate(ci_cross_zero = as.character(ci_cross_zero))%>%
  mutate(ci_cross_zero = ifelse(ci_cross_zero == "0", "yes","no"))%>%
  left_join(param_key2) %>%
  select(param_name2, `2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ2,b2) %>%
  mutate(`2.5%` = round(`2.5%`,3))%>%
  mutate(`97.5%` = round(`97.5%`,3))%>%
  mutate(ci_width = round(ci_width,3))

#Putting together means, credible intervals, and associated stats
EAWP2_summary <- EAWP2_stats %>%
  left_join(EAWP2_quant) %>%
  select(param_name2, Mean, SE_Naive,`2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ2,b2)

EAWP2_summary

write.csv(occ_probs3, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP2_occ_probs.csv')
write.csv(EAWP2_summary, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP2_model_summary.csv')

length(occ_probs2$mean_occ_prob[occ_probs2$mean_occ_prob > 0.95])

length(occ_probs2$mean_occ_prob[occ_probs2$mean_occ_prob > 0.90])

length(occ_probs2$mean_occ_prob[occ_probs2$mean_occ_prob > 0.75])

length(occ_probs2$mean_occ_prob[occ_probs2$mean_occ_prob > 0.50])

length(occ_probs2$mean_occ_prob[occ_probs2$mean_occ_prob > 0.01])

#Retrying original model but without quadratics####

occ_model4 <- jags.model(
  'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_Bayes_Occ3.txt',
  data = list(n = nrow(occ_dat),
              num_tree_spp = occ_dat$num_tree_spp,
              live_basal = occ_dat$live_basal,
              oak_prop = occ_dat$oak_prop,
              can_close_prop = occ_dat$can_clos_prop,
              shrub_dens = occ_dat$shrub_dens,
              prop_forest_1km = occ_dat$prop_forest_1km,
              prop_forest_10km = occ_dat$prop_forest_10km,
              tssr_vis_1 = occ_dat$tssr_vis_1,
              jday_vis_1 = occ_dat$jday_vis_1,
              obs_vis_1 = occ_dat$obs_vis_1,
              tssr_vis_2 = occ_dat$tssr_vis_2,
              jday_vis_2 = occ_dat$jday_vis_2,
              obs_vis_2 = occ_dat$obs_vis_2,
              detect_vis_1 = occ_dat$detect_vis_1,
              detect_vis_2 = occ_dat$detect_vis_2),
  n.chains=3,
  n.adapt=100)

#Making a data frame that is the key for my different parameters####
#This next vector can also be used in my model fit
b3 <- c('b00', 'b01','b02', 'b03','b04', 'b05','b06', 'b07', 'b08','b09', 'b10')
param_name3 <- c("intercept", "tree_spp_richness", "live_basal_area", "oak_prop", "canopy_closure", "shrub_density", "prop_forest_in_1km", "prop_forest_in_10km", "minutes_past_sunrise", 'julian_day','observer')
det_or_occ3 <- c('intercept',rep('occupancy',7),rep('detection',3))

param_key4 <- cbind(b3,param_name3,det_or_occ3) %>%
  as_tibble() %>%
  rename(b = b3, param_name = param_name3, det_or_occ = det_or_occ3)

# sample the MCMC chains 10000 times as burn-in, ignore the samples
update(occ_model4, 10000)

# sample another 10000 times, keep the variables of interest####
EAWP.fit3 <- coda.samples(occ_model4, 
                         c(b3, "true_occ"),
                         10000,
                         thin=1)

#Get summaries of model####

#Mean occupancy probabilities
occ_probs4 <- do.call(rbind.data.frame, EAWP.fit3)[12:512] %>%
  summarize_all(list(mean)) %>%
  matrix(nrow=1) %>%
  t() %>%
  as_tibble() %>%
  rename(mean_occ_prob = V1) %>%
  bind_cols(occ_dat_org) %>%
  select(point_id, mean_occ_prob) %>%
  mutate(mean_occ_prob = unlist(mean_occ_prob))

#Histogram of mean occupancy probs
ggplot(data = occ_probs4, aes(x = mean_occ_prob))+
  geom_histogram(binwidth = 0.1)+
  theme_classic()

#Getting means and standard errors
EAWP_stats3 <- summary(EAWP.fit3)$statistics %>%
  as.data.frame() %>%
  rownames_to_column(var = "b")%>%
  as_tibble()%>%
  filter(str_detect(b, "^b")) %>%
  select(-`Time-series SE`, -SD)%>%
  rename(SE_Naive = `Naive SE`) %>%
  rowwise() %>%
  mutate(Mean = round(Mean, 3))%>%
  mutate(SE_Naive = round(SE_Naive,3)) %>%
  left_join(param_key4) %>%
  select(param_name, Mean, SE_Naive, det_or_occ, b) %>%
  mutate(abs_mean = abs(Mean)) %>%
  arrange(desc(abs_mean)) %>%
  select(-abs_mean)

#Getting credible intervals, determining their width, and seeing if they cross zero
EAWP_quant3 <- summary(EAWP.fit3)$quantiles %>%
  as.data.frame() %>%
  rownames_to_column(var = "b")%>%
  as_tibble() %>%
  filter(str_detect(b, "^b")) %>%
  select(b, `2.5%`, `97.5%`) %>%
  rowwise() %>%
  mutate(ci_width = `97.5%` - `2.5%`) %>%
  mutate(low_pos =ifelse(`2.5%` > 0, 1, 0))%>%
  mutate(up_pos = ifelse(`97.5%` > 0, 1, 0))%>%
  mutate(ci_cross_zero = ifelse(low_pos == up_pos,1,0))%>%
  arrange(desc(ci_cross_zero), ci_width)%>%
  mutate(ci_cross_zero = as.character(ci_cross_zero))%>%
  mutate(ci_cross_zero = ifelse(ci_cross_zero == "0", "yes","no"))%>%
  left_join(param_key4) %>%
  select(param_name, `2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ,b) %>%
  mutate(`2.5%` = round(`2.5%`,3))%>%
  mutate(`97.5%` = round(`97.5%`,3))%>%
  mutate(ci_width = round(ci_width,3))

#Putting together means, credible intervals, and associated stats
EAWP_summary3 <- EAWP_stats3 %>%
  left_join(EAWP_quant3) %>%
  select(param_name, Mean, SE_Naive,`2.5%`, `97.5%`, ci_width, ci_cross_zero, det_or_occ,b)

EAWP_summary3

write.csv(occ_probs4, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_occ_probs_3.csv')
write.csv(EAWP_summary3, 'Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_model_summary_3.csv')


occ_probs <- read_csv('Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP_occ_probs.csv')
occ_probs3 <- read_csv('Data Exploration and Analysis/EAWP Bayesian Occupancy/EAWP2_occ_probs.csv')


length(occ_probs4$mean_occ_prob[occ_probs4$mean_occ_prob > 0.95])
length(occ_probs4$mean_occ_prob[occ_probs4$mean_occ_prob > 0.90])
length(occ_probs4$mean_occ_prob[occ_probs4$mean_occ_prob > 0.75])
length(occ_probs4$mean_occ_prob[occ_probs4$mean_occ_prob > 0.50])
length(occ_probs4$mean_occ_prob[occ_probs4$mean_occ_prob > 0.10])
length(occ_probs4$mean_occ_prob[occ_probs4$mean_occ_prob > 0.01])

length(occ_probs$mean_occ_prob[occ_probs$mean_occ_prob > 0.99])
length(occ_probs$mean_occ_prob[occ_probs$mean_occ_prob > 0.95])
length(occ_probs$mean_occ_prob[occ_probs$mean_occ_prob > 0.90])
length(occ_probs$mean_occ_prob[occ_probs$mean_occ_prob > 0.50])
length(occ_probs$mean_occ_prob[occ_probs$mean_occ_prob > 0.10])
length(occ_probs$mean_occ_prob[occ_probs$mean_occ_prob > 0.05])
length(occ_probs$mean_occ_prob[occ_probs$mean_occ_prob > 0.01])



length(occ_probs3$mean_occ_prob[occ_probs3$mean_occ_prob > 0.95])
length(occ_probs3$mean_occ_prob[occ_probs3$mean_occ_prob > 0.90])
length(occ_probs3$mean_occ_prob[occ_probs3$mean_occ_prob > 0.75])
length(occ_probs3$mean_occ_prob[occ_probs3$mean_occ_prob > 0.50])
length(occ_probs3$mean_occ_prob[occ_probs3$mean_occ_prob > 0.10])
length(occ_probs3$mean_occ_prob[occ_probs3$mean_occ_prob > 0.01])

