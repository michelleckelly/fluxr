library(rstan)
library(devtools)
library(tidyverse)
library(roxygen2)

# Build package
roxygen2::roxygenize()

# Compile stan model
modl <- stanc(file = "inst/stan/o2_twostation.stan")
str(modl)

## Note: for cleaning new datasets: will need to import convert_UTC_to_solartime function from streammetabolizer


# Load internal dataset ########################################################
# Already has equal time breaks, travel time calculated, no NAs
data(diel)

# Pull just the O2 data 
diel <- 
  diel %>%
  select(solar.time, location, DO.obs, DO.sat, temp.water, light.calc,
         depth, discharge, tt, pressure)

# Prep dataset for modeling ####################################################

# Grab number of dates in time series (d)
d <- length(unique(date(diel$solar.time)))
# Convert K600 values to kO2 values?

# Break into dataframes by location
up <- diel[diel$location == "upstream",]
down <- diel[diel$location == "downstream",]

# If we were dealing with multiple days of data, I'd break the dataset up with a for loop here, one loop to model each day
# (for now, select only data from day 2):
d <- 1

# Average daily depth
z_m <- mean(diel$depth)

# Average daily travel time in minutes
tt_min <- mean(diel$tt)


# Seperate datasets by location
up <- up[date(up$solar.time) == ymd("2024-09-19"),]
down <- down[date(down$solar.time) == ymd("2024-09-19"),]

# Data dimensions
#
# Grab timestep (time difference between measurements)
timestep <- as.numeric(up$solar.time[2] - up$solar.time[1]) # min
#
# Calculate n24, max number of observations in 24 hours
n24 <- round(1440 / timestep)
#
# Grab number of observations per date (n)
n <- min(c(length(up$solar.time), length(down$solar.time)))

# Calculate lag using travel time and timestep
lag <- round(tt_min / timestep)

# Calculate light fraction
# Add as column? or vector?
i <- 50
lightstep_up <- sum(up$light.calc[i:(i+lag)])
lighttotal_up <- sum(up$light.calc)
lightfrac_up <- lightstep_up / lighttotal_up

lightstep_down <- sum(down$light.calc[i:(i+lag)])
lighttotal_down <- sum(down$light.calc)
lightfrac_down <- lightstep_down / lighttotal_down

# Avg lightfrac for time window
mean(c(lightfrac_up, lightfrac_down))

# Sample stan model (see runstan_bayes & prepdata_bayes)
#rstan::sampling()



sum(down$light.calc)


# Separate data into vectors
DO.obs.up <- up$DO.obs
DO.obs.down <- down$DO.obs
DO.sat.up <- up$DO.sat
DO.sat.down <- down$DO.sat
# light_mult_GPP?
# const_mult_ER?
# KO2_conv
























# Using travel time and timestep, calculate lag

# Separate into , & only interested in O2






# Pull just the O2 data
