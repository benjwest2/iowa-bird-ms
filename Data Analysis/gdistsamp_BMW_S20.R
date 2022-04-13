
#Load functions from the other file. You'll need tidyverse and unmarked installed to do this.
source("Lunchinators/Abrev gdistsamp functions.R")

#Load data####
#Birds and covariates that can change between visits
birds <- read_csv("Lunchinators/ACFL_data_2018_2019.csv")
#Static (habitat) covariates. These data are center scaled
hab <- read_csv("Lunchinators/veg_scale_2018_2019.csv")


#Look at current data structure
head(birds)
str(birds)

head(hab)
str(hab) #There are some weird "r" column labels in here; don't worry about them

#Get data into a format that unmarked can read. I wrote a custom function that helps do that.
ACFL_data <- make_gdistsamp_df(sp = "ACFL",
                               static_covs = hab,
                               birds_and_dyn_covs = birds, 
                               dcovnames = c("jday_r", "jday_rsq","tssr_r","tssr_rsq","wind_r","obs"),
                               dbreaks = c(0,25,50,75,100))

#Look at current data structure
str(ACFL_data)

#There are four primary pieces of this sort of hierachical distance sample model:####

#1) Detection (aka key) function (What function best models how detection declines with distance?)

#2) Detectability (Given that a bird is available to be detected, does an observer detect it?)
#   -This is assessed using distance bins (detection declines with distance); can also include covariates.

#3) Availability (Given that a bird is using a site, what is the probability that it will make itself available to be detected?)
#   -This assessed using multiple visits, with covariates that differ between visits

#4) Abundance (What affects the number of birds at a site?)
#   -This is assessed using site (i.e., habitat) covariates

#######################
#Detection function####
#######################

#First, picking a detection (aka key) function using an intercept-only model. I wrote a function to do this.
#This takes about 5 minutes on my laptop
keyfuns <- keyfun_choose(ACFL_data, site_covs = 1, det_covs = 1, av_covs = 1, funs = c("halfnorm", "hazard", "exp"))

str(keyfuns)

#best key function. In this case, it is a half normal
keyfuns[[1]]

#AIC table
keyfuns[[2]]

#############################
#Detectability covariates####
#############################

#We'll just run all combinations; it's just four models. Observer ("obs") is a factor with three levels. Wind speed ("wind") has a 0 for minimum wind speed and 1 for maximum observed wind speed and is continuous between.

#Null model
null_mod <- gdistsamp(data = ACFL_data,~1, ~1, ~1, keyfun = "halfnorm", output = "density", unitsOut = "ha")

#Observer only for detection
det_obs <- gdistsamp(data = ACFL_data, ~1, ~1, ~obs, keyfun = "halfnorm", output = "density", unitsOut = "ha")

#Wind only for detection
det_wind <- gdistsamp(data = ACFL_data, ~1, ~1, ~wind, keyfun = "halfnorm", output = "density", unitsOut = "ha")

#Wind + observer for detection
det_all <- gdistsamp(data = ACFL_data, ~1, ~1, ~wind + obs, keyfun = "halfnorm", output = "density", unitsOut = "ha")

#Check detectability models with AIC
det_AIC <- modSel(fitList(intercept = null_mod, obs = det_obs, wind= det_wind, "wind + obs" = det_all))

#Display AIC table for detectability
det_AIC


#getting best detectabilty model and formatting it as a one-sided formula
best_det <- paste("~",det_AIC@Full$model[1], sep = "")


############################
#Availability covariates####
############################

#For the sake of time, I'm only going to run 3 different models, but I have run all combinations previously. tssr is time since sunrise, and jday is Julian day. _sq indicates the variable is squared. These values are standardized such that the maximum value is 1 and the minimum is 0.

#Null model
av_null <- det_obs

#Just tssr and jday (no quadratics)
av_no_quad <- gdistsamp(data = ACFL_data,~1, ~tssr + jday, best_det, keyfun = "halfnorm", output = "density", unitsOut = "ha")

#tssr and jday with quadratics
av_all <- gdistsamp(data = ACFL_data, ~1, ~jday + jday_sq + tssr + tssr_sq, best_det, keyfun = "halfnorm", output = "density", unitsOut = "ha")

#Check availability models with AIC
av_AIC <- modSel(fitList("intercept only" = av_null, "tssr + jday" = av_no_quad, "tssr + tssr_sq + jday + jday_sq" = av_all))

#Display AIC table for availability
av_AIC

#getting best availabilty model and formatting it as a one-sided formula
best_av <- paste("~",av_AIC@Full$model[1], sep = "")

#########################
#Abundance covariates####
#########################

#Now running a few models with different abundance (/density) covariates, using the best combination of detectability and availability covariates.
#I go over how I actually do model selection for this in my presentation slides, as well as what the different covariates are 

#No abundance covariates
abund_null <- av_all

#A priori hypothesized model
abund_pred <- gdistsamp(data = ACFL_data,
                        ~nearest_patch_size + fprop_1km + fprop_10km + 
                          dist_edge + 
                          live_basal,
                        best_av, 
                        best_det, 
                        keyfun = "halfnorm", 
                        output = "density", 
                        unitsOut = "ha")

#Top model (determined previously)
abund_top <- gdistsamp(data = ACFL_data,
                       ~nearest_patch_size + fprop_10km + 
                         dist_edge + 
                         spp_rich + live_basal + oak_prop + 
                         grass_prop + litter_prop + shrub_dens, 
                       best_av, 
                       best_det, 
                       keyfun = "halfnorm", 
                       output = "density", 
                       unitsOut = "ha")
#Full model
abund_full <- gdistsamp(data = ACFL_data,
                       ~year +
                         nearest_patch_size + fprop_1km + fprop_10km + 
                         dist_edge + 
                         spp_rich + live_basal + dead_basal + oak_prop + 
                         grass_prop + litter_prop + green_prop + shrub_dens, 
                       best_av, 
                       best_det, 
                       keyfun = "halfnorm", 
                       output = "density", 
                       unitsOut = "ha")

#AIC for abundance models
abund_AIC <- modSel(fitList("intercept only" = abund_null, "top model" = abund_top, "predicted model" = abund_pred, "full model" = abund_full))

#Display AIC table for abundance models
abund_AIC

####################################
#Summarize the top model results####
####################################

#I have more veg points than usable bird points; fixing that
hab_trim <- hab %>%
  filter(point_yr %in% birds$point_yr)

#Now for an actual summary. This uses another custom function that wraps some unmarked base functions
top_mod_summary <- summarize_gds(abund_top, hab_trim, point_cols = c("bca","point_yr","year"), remove_suffix = T, count_radius = 100)

#Look at summary structure
str(top_mod_summary)

#top_mod_summary$point_ests is estimates for each point_year combination

#Summaries of coefficients by process####

#Detectability
top_mod_summary$coeff_summary%>%
  filter(process == "detectability")

#Availability
top_mod_summary$coeff_summary%>%
  filter(process == "availability")

#Abundance
top_mod_summary$coeff_summary%>%
  filter(process == "abundance")

#Rough density estimates at a site-year level

counts_by_site <- top_mod_summary$point_ests %>%
  mutate(indicator = 1)%>%
  group_by(bca, year)%>%
  summarize(numbirds = sum(mean_count), sum_low95ci = sum(count_ci_95_low), sum_up95ci = sum(count_ci_95_up), num_points = sum(indicator))

#Note that the "credible intervals" here are just sums of the point scale ci's; the lower sum of ci's is eqaul to or slightly larger than the observed count
counts_by_site

#Convert counts to densities. Same caveats about the "ci's"
counts_by_site %>%
  #this is area in hectares
  mutate(point_area = pi*100^2*0.0001)%>%
  mutate(survey_area = num_points*point_area)%>%
  mutate(mean_dens = numbirds/survey_area)%>%
  mutate(sum_low95ci_dens = sum_low95ci/survey_area)%>%
  mutate(sum_up95ci_dens = sum_up95ci/survey_area)%>%
  select(bca, year, mean_dens, sum_low95ci_dens, sum_up95ci_dens)


