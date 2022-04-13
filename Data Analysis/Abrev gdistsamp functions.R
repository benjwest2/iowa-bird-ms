#A slew of functions to run gdistamp models in unmarked. (This script isn't
#standalone and is meant to be called y "gdistsamp_BMW_S20.R")

#Libraries####
library(unmarked)
library(tidyverse)

#A function to take distance sampling bird point count data and habitat data and turn it into an unmarked data frame usable by gdistsamp.####

#First argument is species ('sp'), which is a string corresponding to a banding code (e.g., "EAWP"). This code must exist in 'birds_and_dyn_covs'.

#Second argument ('static_covs') is a data frame of site covariates that the is only one values for, e.g., field-collected vegetation data. A column with prefix "point" (e.g., "point_id", "point_yr") is required, and each covariate gets a single column.

#The third argument ('birds_and_dyn_covs') contains bird counts, and covariates that vary between years and/or visits; this data frame must have very specific formatting. All covariates other than observer (or 'obs') must be numeric; if not, code will need to be modified. I will define "n" as the number of visits to site. The first two columns should be 'point_id' and 'species', where species has 4-letter banding codes. The next columns should be nbins*n columns of point count data stratified into nbins number distance bins. The primary ordering needs to be by visit number, then secondarily by distance bin number (e.g. for n=2, nbins=2, bin1_vis1, bin2_vis1, bin1_vis2, bin2_vis2). The next columns must be n replicates of any number of covariates. For each set of n covariate columns, each column in the set must have the same prefix ending in a "_", e.g., "year_8". This prefix must be unique to a set. The columns of covariates in a set must be ordered chronologically by visit. The prefixes of one of these sets must start with "obs" and correspond with observer, which will be turned into a factor.

#All "point" column values in 'birds_and_dyn_covs' must have a match in 'static_covs', and the column names must be identical. 

#Fourth argument ('dcovnames') is a vector of character strings corresponding to prefixes of covariate column names in 'birds_and_dyn_covs'. All individual strings should omit the final "_" , e.g. c("year", "obs") NOT c("year_", etc.,) and must occur as prefixes of sets of column names in 'birds_and_dyn_covs'. I have a default set as c("jday_r", "jday_rsq","tssr_r","tssr_rsq","wind_r","obs","fprop_1km","fprop_10km").

#Fifth argument: dbreaks is a vector of the radii of the different count bands in meters, with the additon of 0 as the first value (defaults to c(0,25,50,75,100)). 

make_gdistsamp_df <- function(sp, static_covs, birds_and_dyn_covs, dcovnames = c("jday_r", "jday_rsq","tssr_r","tssr_rsq","wind_r","obs"), dbreaks = c(0,25,50,75,100)){
  
  require(unmarked)
  require(tidyverse)

  #First figuring out what the column with the "point" prefix is called
  
  point_colname <- static_covs %>%
    #Filter to species of interest
    select(starts_with("point"))%>%
    names()%>%
    unlist()
  
  #Upfront processing to birds_and_dyn_covs that will used for both the bird and covariate data
  bird_yr_cov2 <- birds_and_dyn_covs %>%
    #Filter to species of interest
    filter(species == sp)%>%
    #Making sure ordering of points is consistent
    arrange(!!sym(point_colname)) 
  
  #Getting vector of point/point-year ids that have bird data. Important because I have habitat data from all points, but bird data from some points
  point_vect <- bird_yr_cov2 %>%
    select(!!sym(point_colname)) %>%
    unlist(use.names = F)
  
  #Isolate bird counts to a matrix
  bird_counts <- bird_yr_cov2 %>%
    #Just get point_id and distance binned counts
    select(!!sym(point_colname), starts_with("dist"))%>%
    #Make point_id row names and convert to matrix
    rename(rowname = !!sym(point_colname)) %>%
    column_to_rownames() %>%
    as.matrix()
  
  #Calculate the number of visits to each point
  num_visits = dim(bird_counts)[2]/(length(dbreaks)-1)
  
  #Process static covariates to a vanilla data frame
  sitecovs <- static_covs %>%
    #Paring down static covariates to just points that have bird data
    filter(!!sym(point_colname) %in% point_vect)%>%
    #Making sure ordering of points is consistent
    arrange(!!sym(point_colname))%>%
    #Add row names and convert to vanilla data frame. 
    #The unclass bit makes characters into factors by passing the unclassed data through as.data.frame's defaults
    unclass()%>%
    as.data.frame()%>%
    rename(rowname = !!sym(point_colname)) %>%
    column_to_rownames()
  
  #Need to make year a factor if it is present as a static covariate; i.e., when each row is a point-year
  
  if(TRUE %in% grepl("year", names(sitecovs))){
    
    #Getting year column name. Should be "year" but playing it safe
    year_colname <- static_covs %>%
      #Filter to species of interest
      select(starts_with("year"))%>%
      names()%>%
      unlist()
    
    sitecovs <- sitecovs %>%
      mutate(year = as.factor(!!sym(year_colname)))
  }
    
    
    
  #Process dynamic covariates into a list of data frames
  
  #Doing observer and year stuff before proceeding because they are factors with multiple levels, but not all columns contain all levels 
  
  #Starting with a conditional for observer so I don't get errors if there isn't an observer column
  if(TRUE %in% grepl("obs", names(bird_yr_cov2))){
  #Getting factor levels for observer

  obs_levels <- bird_yr_cov2 %>%
    #Just get observers
    select(starts_with("obs"))%>%
    #Data frame to vector
    unlist(use.names = F) %>%
    #Get rid of duplicates
    unique()
    
  obs_df <- bird_yr_cov2 %>%
    #Get just point_id and observers
    select(!!sym(point_colname), starts_with("obs"))%>%
    #Convert to vanilla data frame
    as.data.frame()%>%
    rename(rowname = !!sym(point_colname)) %>%
    column_to_rownames()
  
  #Setting consistent factor levels across all observer columns. Couldn't figure out how to get this in the pipe without issues.
  obs_df[] <- lapply(obs_df[], factor, 
                        levels=c(obs_levels), 
                        labels = c(obs_levels))
  
  #Remove observer from covariate name list
  dcovnames <- dcovnames %>%
    .[-grep("obs", .)]
  }
  #END OBSERVER CONDITIONAL####
  
  #Year stuff for when year is a dynamic covariate for each point (not point-year)####
  
  #Starting with a conditional for year in bird_yr_cov2 so I don't get errors if there isn't a year column
  if(TRUE %in% grepl("year", names(bird_yr_cov2))& "year" %in% dcovnames){
    #Getting factor levels for observer
    
    year_levels <- bird_yr_cov2 %>%
      #Just get observers
      select(starts_with("year"))%>%
      #Data frame to vector
      unlist(use.names = F) %>%
      #Get rid of duplicates
      unique()%>%
      as.character()
    
    year_df <- bird_yr_cov2 %>%
      #Get just point_id and observers
      select(!!sym(point_colname), starts_with("year"))%>%
      #Convert to vanilla data frame
      as.data.frame()%>%
      mutate_all(as.character())%>%
      rename(rowname = !!sym(point_colname)) %>%
      column_to_rownames()
    
    #Setting consistent factor levels across all observer columns. Couldn't figure out how to get this in the pipe without issues.
    year_df[] <- lapply(year_df[], factor, 
                       levels=c(year_levels), 
                       labels = c(year_levels))
    
    #Remove year from dynamic covariate name list
    dcovnames_yr <- dcovnames
    dcovnames <- dcovnames %>%
      .[-grep("year", .)]
  }else{
    dcovnames_yr <- NULL
  }
  
  #END DYNAMIC YEAR CONDITIONAL####
  
  #Making a cleaner version of dcovnames
  dcovnames_clean <- dcovnames %>%
    str_replace("_r","") %>%
    str_replace("sq","_sq")
    
  #Making an empty list of data frames to fill with a for loop
  dcov_list <- vector(mode = "list", length = length(dcovnames))
  
  #For loop to fill a list with data frames. Each data frame is a dynamic covariate that varies between visits or years
  for (i in 1:length(dcovnames)){
    cov <- dcovnames_clean[i] %>%
      #Need to add the "_" because otherwise jday and tssr pick up jday_sq and tssr_sq
      assign(x = ., value = select(bird_yr_cov2,starts_with(paste((dcovnames[i]),"_", sep=""))))%>%
      #Needs to be a vanilla data frame
      as.matrix()
    dcov_list[[i]] <- cov
  }
  
  #Assigning names to the dynamic covariate list of data frames
  dcov_list_final <-  dcov_list
  names(dcov_list_final) <- dcovnames_clean
  
  #Adding in observers if they were present
  if(TRUE %in% grepl("obs", names(bird_yr_cov2))){
  dcov_list_final <- append(dcov_list_final,list(obs_df))
  names(dcov_list_final)[length(dcov_list_final)] <- "obs"
  }
  
  #Adding years if they were present as a dynamic covariate
  if(TRUE %in% grepl("year", names(bird_yr_cov2))&"year" %in% dcovnames_yr){
    dcov_list_final <- append(dcov_list_final,list(year_df))
    names(dcov_list_final)[length(dcov_list_final)] <- "year"
  }
  
  
  umf <- unmarkedFrameGDS(y = bird_counts, siteCovs = sitecovs, numPrimary = num_visits, yearlySiteCovs = dcov_list_final, dist.breaks = c(0,25,50,75,100), survey = "point", unitsIn = "m")
  
  return(umf)
}


#A function to make a character vector into a character than can be used as the right hand side of a regression equation. Example: c("obs","wind") becomes "~obs + wind". Setting tilde = FALSE will not put the tilde in front of the formula; default is with tilde.

vect_to_reg_eq <- function(char_vect, tilde = TRUE){
  
  reg_eq <- char_vect %>%
  paste(collapse = " + ")%>%
  ifelse(tilde == TRUE, paste("~",.,sep = ""),.)
  
  return(reg_eq)
}

#A function to make an AIC table given a character vector or data frame column that corresponds to names of 'unmarkedFitGDS' objects ("models"), and an optional data frame ("more_info_df") that has model names (name must be "mod_names") as column 1 and any additional info about the models as other columns. The output is a tibble, and it is arranging in terms of ascending AIC, with the lowest AIC model as row 1.


AIC_table_GDS <- function(models, more_info_df = "NONE"){
  
  AICs <- c()
  #Make sure "models" is vector
  models <- unlist(models, use.names = F)
  
  #Make a vector of AIC tables
  for(i in 1:length(models)){
    model <- get(models[i])
    AIC_ <- model@AIC
    AICs[i] <- AIC_
  }
  
  #Make an AIC table
    AIC_table <- tibble(mod_name = models, AIC = AICs) %>%
    mutate(AIC = as.numeric(AIC)) %>%
    mutate(dAIC = AIC - min(AIC)) %>%
    arrange(AIC)%>%
    mutate(AIC = as.character(round(AIC,2)))%>%
    mutate(dAIC = as.character(round(dAIC,2)))
    
  #Add model formulas to table if specified
    if(is.data.frame(more_info_df)){
      AIC_table <- AIC_table %>%
        left_join(more_info_df)
    }
    
    
    return(AIC_table)
  
}

#A function to choose the top detection key function (half normal, negative exponential, hazard, uniform).It returns a list with two elements: The first is the name of the best fitting function, the second is the an AIC table of the different function fits
#It requires:
#1) unmarked data frame ("umf")
#2) a vector of site covariates ("site_covs"; "year" can be included if it is the same for all visits in a row; defaults to all site covariates in the unmarked data frame). All names in this vector must occur in umf's siteCovs section
#3) a vector of detectability covariates ("det_covs",defaults to c("obs","wind"). All names in this vector must occur in umf's yearlySiteCovs section
#4)a vector of covariate names that affect availability ("av_covs"; defaults to c("tssr","tssr_sq","jday","jday_sq")). All names in this vector must occur in umf's yearlySiteCovs section
#5) You can optionally choose fewer key functions using 'funs'
#REQUIRES CUSTOM FUNCTION vect_to_reg_eq TO TURN vectors into the right hand side of a regression equation

keyfun_choose <- function(umf, site_covs = "ALL", det_covs = c("obs","wind"), av_covs = c("tssr","tssr_sq","jday","jday_sq"), funs = c("halfnorm", "hazard", "exp", "uniform")){
  
  require(unmarked)
  require(tidyverse)

  #get covariates vectors into regression equation format
  av <- vect_to_reg_eq(av_covs)
  det <- vect_to_reg_eq(det_covs)
  
  if(site_covs == "ALL"){
    site <- umf@siteCovs%>%
      names()%>%
      vect_to_reg_eq()
  } else {
    site <- vect_to_reg_eq(site_covs)
  }
  
  
  #K is number of sites; easiest way to get that is to pull it from umf
  #Run each model, and assign each model to a variable with the name of the function  
  for(i in 1:length(funs)){
    assign(funs[i], gdistsamp(site, phiformula = av, pformula = det, umf, output="density", keyfun = funs[i], K=length(umf@tlength)), env = .GlobalEnv)
    print(paste("Function model",i,"of",length(funs),"complete:",funs[i],"at", Sys.time()))
  }
  
  #Make an AIC table
  AIC_table <- AIC_table_GDS(funs)
  
  #Get name of best model
  best_mod_name <- AIC_table[1,1] %>%
    unlist(use.names = F)
  
  return(list(best_key_fun = best_mod_name, key_fun_AIC = AIC_table))
}


#A function to summarize the results of a gdistsamp analysis. It returns a list of tibbles: the first tibble ("point_ests") is point-scale estimates including mean density in hectares, the mean empirical Bayes estimate for number birds within a point's count radius, and credible intervals for the number of birds within a point's count radius (95% is default). The second tibble ("coeff_summary") includes estimates of coefficient values and 95% confidence intervals for all covariates on lambda (abundance), phi (availability), and p (detectability), as well as intercepts for each parameter.

#The first argument (gds_mod) is an output of a gdistsamp model (class = 'unmarkedFitGDS). The second argument (point_df) is a data frame containing a single column whose label has the prefix "point" with exactly the same points that were used in the gdistsamp model. If you used gdist_assist function "make_gdistsamp_df," this will be the same data frame that you used for static_covs. The third argument ("point_cols") is a character vector of column names within point_data to include in the point estimate output, written in the order they will appear in the output. The third argument ("remove_suffix", default = FALSE) will, if set to TRUE, remove a suffix from point IDs within "point_data" , given there is a single separator (e.g., "_", ".", " ") between the point ID and the suffix. An example would be pointID_2020, where 2020 is the suffix to remove. The fourth argument ("point_ci") is a number between 0 and 100 specifying the desired % credible interval cutoff for point_ests; default is 95. The fourth argument ("coeff_ci" ) is a number between 0 and 100 specifying the desired % confidence interval cutoff for coeff_summary; default is 95. The fifth argument ("count_radius) is the truncated point count radius in meters.


#WARNING: Make sure point names are sorted alphabetically (A-Z) in your unmarked data frame before running your gdistsamp model. The gdist_assist function "make_gdistsamp_df" does this automatically.

summarize_gds <- function(gds_mod, point_df, point_cols, remove_suffix = FALSE, point_ci = 95, count_radius){
  
  #Estimate point scale densities####
  
  #Get point column name
  point_colname <- point_df %>%
    select(starts_with("point"))%>%
    names()%>%
    unlist()
  
  if (!(point_colname %in% point_cols)){ 
    stop("point column name in 'point_df' does not match point column name in 'point_cols'")
  }
  
  
  #Get point year combinations to assign to estimates
  point_data_trim <- point_df %>%
    select(point_cols)%>%
    arrange(!!sym(point_colname))
  
  if(remove_suffix==T){        
    point_data_trim <-  point_data_trim %>%
      separate(!!sym(point_colname), into = c("point","junk"))%>%
      select(-junk)
    
    point_cols <- replace(point_cols, point_cols == point_colname, "point")
    
  }
  
  #Get empirical Bayes estimates of # birds for each point####
  pt_ranef <- ranef(gds_mod)
  
  point_ests <- pt_ranef %>%
    #Note that these are upper and lower credible intervals for COUNTS, not density
    confint(level = (point_ci/100))%>%
    as_tibble()%>%
    rename(count_ci_95_low = `2.5%`, count_ci_95_up = `97.5%`)%>%
    #bup = best unbiased estimator. Used to extract mean count at a point.
    mutate(mean_count = bup(pt_ranef, stat = "mean"))%>%
    mutate(mean_dens = mean_count/(pi*count_radius^2)*10000)%>%
    bind_cols(point_data_trim)%>%
    select(one_of(point_cols), mean_dens, mean_count, count_ci_95_low, count_ci_95_up)
  
  #Get parameter estimates and confidence intervals####
  
  #parameter names
  parameters <- names(gds_mod@estimates@estimates)
  processes <- c("abundance", "availability", "detectability")
  
  #blank tibble 
  coeff_summary <- matrix(ncol=6, nrow=0, dimnames= list(NULL,c("process","parameter","covariate","mean_coeff", "ci_95_low", "ci_95_up"))) %>%
    as_tibble()%>%
    mutate_at(c("process","parameter","covariate"), as.character)%>%
    mutate_at(c("mean_coeff", "ci_95_low", "ci_95_up"), as.numeric)
  
  
  for(i in 1:length(parameters)){
    coeff_summary <- gds_mod %>%
      confint(type = parameters[i]) %>%
      as.data.frame()%>%
      rownames_to_column(var = "coeff")%>%
      separate(coeff, into = c("parameter","covariate"), sep = "\\(")%>%
      mutate(covariate = gsub("\\)", "", covariate))%>%
      mutate(covariate = gsub("Int", "(intercept)", covariate))%>%
      mutate(process = processes[i])%>%
      rename(ci_95_low = `0.025`,  ci_95_up =`0.975`)%>%
      as_tibble()%>%
      mutate(mean_coeff = gds_mod@estimates@estimates[[parameters[i]]]@estimates)%>%
      select(process,parameter,covariate,mean_coeff, ci_95_low, ci_95_up, process)%>%
      bind_rows(coeff_summary)
  }
  
  return(list(point_ests = point_ests, coeff_summary = coeff_summary))
  
}



