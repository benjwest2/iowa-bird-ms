#Script written in 2022 to help explain occupancy modeling in unmarked to a
#friend. Uses my current (at the time of writing) approach to data wrangling,
#which I hope is more elegant than some of the other scripts in this project!

library(tidyverse)
library(unmarked)
library(skimr)
 
#Sampling scheme here was two visits per year to each point with distance
#sampling. We're just gonna do a single species, single season occupancy model.

#Start by reading in data####

#Birds and covariates that can change between visits
birds <- read_csv("Input Data/ACFL_data_2018_2019.csv")

#Static (habitat) covariates. These data are center scaled, which means
#they were transformed to all have a mean of 0 and standard deviation of 1 
#(don't worry about it)
hab <- read_csv("Input Data/veg_scale_2018_2019.csv")


skim(birds)
#Lowdown on columns:
#point_yr is unique ID
#the "dist" columns are distance bins (which need to be compressed)
#the "_1" and "_2" suffixes denote which visit the variable corresponds to

skim(hab)
#point_yr is still unique ID
#don't worry about what the other habitat variables are

#Data cleaning####

#Gotta do a little processing... these are distance sampling data. We need
#observed vs. not observed. Also we only want one year of data to satisfy
#closure assumption of single season occupancy modeling. Also unmarked data 
#frames are a pain

#First, what does unmarked need from us to make its flavor of data frame?
?unmarkedFrameOccu

#We need 
#1) y (observations in matrix format)
#2) siteCovs, which are generally habitat covariates and DO NOT CHANGE 
#   between visits
#3) obsCovs, which are covariates that affect whether or not an animal is 
#   observed


#Let's start with y
y_occ <- 
  birds%>%
  #Let's only do 2019
  #("grepl" matches text, filter pulls relevant rows)
  filter(grepl("2019",point_yr))%>%
  #Arranging by point_yr JUST TO BE SAFE; everything HAS to be in the same order
  #in each part of the unmarked data frame
  arrange(point_yr)%>%
  #we only have one species, so don't need to worry about that
  #we only need the "dist" columns (select pulls relevant columns)
  select(starts_with("dist"))%>%
  #This next line finds columns ending in "_1" and sums each row and creates
  #new column "vis_1" from those sums
  #"mutate" is an easy way to create a new column
  mutate(vis_1 = rowSums(select(.,ends_with("_1"))))%>%
  #Same deal, but with second visit
  mutate(vis_2 = rowSums(select(.,ends_with("_2"))))%>%
  #Ditch distance columns
  select(starts_with("vis"))%>%
  #Great, so distance bins are collapse into a single... count
  #We need presence absence
  #Since we're doing the same thing to everything and don't need new columns,
  #we can mutate across everything using the same function
  #we want each row treated separately, so we'll tell R to do it row-wise
  rowwise()%>%
  mutate(across(everything(),
                #This function replaces anything bigger than 1 with a 1
                function(x){
                  min(c(x,1))
                }
                ))%>%
  #Cool, this is presence/absence but it's not a matrix
  as.matrix()
  #Now it is 
  #(everything in a matrix has to be the same type of data; you can't
  #mix text and numbers like you can in a data frame)

#Just basking in the glory of The Matrix
head(y_occ)

#Now on to siteCovs... this is pretty straight forward
hab_occ <-
hab%>%
  #Just 2019 again
  filter(grepl("2019",point_yr))%>%
  #There are more sites with habitat surveys than bird surveys... fixing that
  filter(point_yr %in% birds$point_yr)%>%
  #arrangeby unique ID
  arrange(point_yr)%>%
  #I'm gonna get rid of a few things we aren't gonna use
  select(-c(point_yr, year, nearest_patch_size, 
            fprop_10km, ends_with("sq"), dead_basal))%>%
  #bca (Bird Conservation Area) needs to be factor, which is text stored as
  #number
  mutate(across(bca, as.factor))%>%
  #make sure we're dealing with a plain data frame
  as.data.frame()
  
head(hab_occ)

#Aight so this one is fun... obsCovs needs to be a list of data frames

#first need to remember what the covariates are called

names(birds)
#Cool, so we only need things right of "year"
#the "_r_" just means rescaled
#the "_rsq_" mean rescaled and squared
#We're not gonna worry about the squared stuff

obs_occ <-
#Let's make a named list!
list(
  #pattern is:
  #name = value
  #I'm just using column prefixes as my values
  #jday is Julian day
  jday = "jday_r_",
  #tssr is time since sunrise
  tssr = "tssr_r_",
  #wind is wind speed
  wind = "wind_r_"
  #we'll do observer later
)%>%
  #Aight, so we're gonna use "lapply," which means "do this for everything in
  #the list
  lapply(
    #lapply needs a function
    function(x){
      #gonna use "birds"
      birds%>%
        #Just 2019 again
        filter(grepl("2019",point_yr))%>%
        #arrange for safety
        arrange(point_yr)%>%
        #just selecting the covariate of interest, which is specified by
        #the list value
        select(starts_with(x))%>%
        #and make sure it's a plain data frame
        as.data.frame()
    }
  )
  
#So what did we just do?

#Check to see if it has names
names(obs_occ) #cool

#Check the structure
str(obs_occ) #list of three data frames, great

#We can call each data frame by name, too, using double brackets!
head(obs_occ[["tssr"]])

#This also works, but is a little dirtier
head(obs_occ$jday)

#Now for observer... this next code block is 
#just to see who we have as observers
{birds%>%
  #Just 2019 again
  filter(grepl("2019",point_yr))%>%
  #just selecting observer
  select(starts_with("obs"))%>%
  #turn into a vector
  unlist()%>%
  #what's in the vector?
  table()
}

#OK, we have "CL" and "JAN." We need to make that a factor, but each column
#needs to have the same factor "levels" (i.e., same text needs to correspond
#to the same number)

#(we're gonna add the factored data frame to the "obs_occ" list as an item
#named "obs")

obs_occ[["obs"]] <-
birds%>%
  #Just 2019 again
  filter(grepl("2019",point_yr))%>%
  #just selecting observer
  select(starts_with("obs"))%>%
  #Now make everything a factor, and force the levels to be the same
  mutate(across(everything(),factor, levels = c("CL","JAN")))%>%
  #and plain data frame time
  as.data.frame()

#observer is in there now!
names(obs_occ)

str(obs_occ)

#Hard part is over; now for assembly!

occ_df <- unmarkedFrameOccu(
  y = y_occ,
  siteCovs = hab_occ,
  obsCovs = obs_occ
)

#Based on that warning, apparently unmarked did the factor conversion again... 
#oh well

#Modeling~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#Two ways of doing the simplest possible occupancy model ("null" model)####
?occu

#Standard
occu(~1~1, data = occ_df)

#Using text then converting to formula (this will become very useful very soon)
occu(as.formula("~1~1"), data = occ_df)

#Stuff after first squiggle is detection/observation, stuff after second
#squiggle is occupancy (based on habitat, is the critter there or not?)

#Detection covariates####

#Let's find the best detection covariates, using ALL of the habitat covariates
#as the baseline for occupancy (you can also do null (~1), by I do it this way)

#first, define the full model for occupancy (I have a dirty trick here)

(full_occ_form <-
hab_occ%>%
  #take names of covariates
  names()%>%
  #paste them together with a " + " separator
  paste(collapse = " + ")%>%
  #put the squiggle up front
  paste0("~",.))

#Make a named list of models

detect_mod_names <-
#We'll keep it pretty simple, null model, full model, and single covariates
  list(null = "~1",
     full = "~jday + tssr + wind + obs",
     jday = "~jday",
     tssr = "~tssr",
     wind = "~wind",
     obs = "~obs")


detect_mod_list <-
detect_mod_names %>%
  #Back to lapply so turn everything into a model
  lapply(function(x){
    occu(as.formula(
      #using list items as first model, then the full model
      paste0(x,full_occ_form)
      ), data = occ_df)
  })
  

#There's some warning here for me, ignore it

#Time to rank some models!

#We'll use AIC; you'll learn more about this if you haven't already,
#but all we need to know right now is that low score = good

#unmarked has a default fancy AIC table, but it's hard to use. This is quick and
#dirty and works just fine for right now

detect_mod_list%>%
  lapply(function(x){
    x@AIC%>%
      data.frame(AIC = .)
  })%>%
  bind_rows(.id = "model")%>%
  #arrange puts low score on top
  arrange(AIC)

#Looks like full model worked best! We'll use that!

#Here's a summary
detect_mod_list[["full"]]

#Doesn't actually look like jday does anything, but closer to sunrise,
#low wind, and having CL (wonder who that is...) as the bird observer gives 
#us our best shot at detecting an Acadian Flycatcher! We *could* remove jday, 
#but really, we just needed more models or better models at the start

#Speaking of models at the start... let's mess with the occupancy covariates

#Occupancy covariates####

#Usually, you always want the the full model and null model in your model set

#Some people will test all possible models, but that's bad form; let's pick
#some models

#1 will be null
#2 will be full
#3: let's do a "landscape model". This will just have 
#   Bird Conservation Area (bca), fprop_1km (proportion of land forested in a 
#   1 km radius), and dist_edge (distance to edge). Idea here is that bird 
#   mostly uses large-scale cues to find a habitat
#4: now a "site vegetation structure" model: live_basal and shrub_dens, 
#   so tree density and shrub density.
#5: "all site qualities" mix of structure and other site qualities, but no
#   landscape stuff. Basically everything but bca, fprop_1km, and dist_edge
#6: Acadian flycatchers are interior forest species, so characteristics interior
#   forest would be large distance to edge, lots of forest on the landscape, few
#   shrubs, and little grass.

#Now the maths


occ_mod_names <-
  #We'll keep it pretty simple, null model, full model, and single covariates
  list(null = "~1",
       full = full_occ_form,
       landscape = "~bca + fprop_1km + dist_edge",
       site_struct = "~live_basal + shrub_dens",
       #this gets a little long; gonna break it up between two lines with paste0
       site_full = paste0("~live_basal + shrub_dens + oak_prop  + ",
                          "grass_prop + green_prop + litter_prop + spp_rich"),
       int_forest = "~fprop_1km + dist_edge + grass_prop + shrub_dens")


occ_mod_list <-
  occ_mod_names %>%
  #Back to lapply so turn everything into a model
  lapply(function(x){
    occu(as.formula(
      #using list items as first model, then the full model
      paste0("~jday + tssr + wind + obs",x)
    ), data = occ_df)
  })
  
#Simple AIC table
occ_mod_list%>%
  lapply(function(x){
    x@AIC%>%
      data.frame(AIC = .)
  })%>%
  bind_rows(.id = "model")%>%
  #arrange puts low score on top
  arrange(AIC)

#So much for the models I picked... full model reigned supreme by a long shot

#Here's a summary
occ_mod_list[["full"]]

#SAND_CREEK BCA is tied up in the intercept, but STEPHENS has fewer sites with
#ACFL, and THOUSAND_ACRES has more sites than either of the other two.

#tree species richness didn't really matter
#more basal area (more trees) = better. Makes sense; some of the points aren't
#forested

#More oak = higher occupancy probablity

#Less grass; makes sense for an interior forest bird

#Green and litter cover don't matter

#fewer shrubs = higher occ prob; makes sense for an interior forest bird

#fprop_1km didn't matter

#distance to edge had a HUGE effect

#One advantage of center scaling is that because every variable is now on the
#scale, you just look at the estimates and say "dist_edge" had a bigger effect
#than "oak_prop"









