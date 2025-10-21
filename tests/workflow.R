library(rstan)
library(devtools)
library(tidyverse)
library(roxygen2)
library(rstan)

# load package
#roxygen2::roxygenize()
load_all()

# Compile stan model
#modl <- stanc(file = "inst/stan/o2_twostation.stan")
#str(modl)

## Note: for cleaning new datasets: will need to import convert_UTC_to_solartime function from streammetabolizer

# Load internal dataset ########################################################
# Already has equal time breaks, travel time calculated, no NAs
data(diel)

# Pull just the O2 data 
data_all <- 
  diel %>%
  select(solar.time, location, DO.obs, DO.sat, temp.water, light.calc,
         depth, discharge, tt, pressure) %>%
  rename(light = light.calc)

# Add date column
data_all$date <- date(data_all$solar.time)

# Set specifications for modeling ##############################################
mod_specs <- 
  set_specs(model_name = "o2_twostation", verbose = TRUE,
            K600_daily_meanlog = 1, K600_daily_sdlog = 0.05,
            n_cores = 8)

# Pass data and specs to prepping function #####################################
data_prepped <- 
  prepdata_bayes_twostation(data = data_all, ply_date = "2024-09-19", 
                            specs = mod_specs)

# Run modeling function ########################################################
mod_out <- 
  runstan_bayes_twostation(data_list = data_prepped, 
                           specs = mod_specs)

# Evaluate model output ########################################################

# Plot modeled O2 vs measured O2
## Function glues together modeled O2 and observed O2 dataframes
return_DO <- function(data_list, model_output){
  
  data.frame(
    date = date(model_output$inst$solar.time),
    solar.time = model_output$inst$solar.time,
    DO.obs.down = data_list$DO_obs_down,
    #DO.sat.down = data_list$DO_sat_down,
    #depth = data_list$depth,
    #temp.water = data_list$temp,
    #light = data_list$light_mult_GPP,
    DO.mod.down = model_output$inst$DO_mod_down_mean,
    DO.mod.down.lower = model_output$inst$DO_mod_down_mean - 
      model_output$inst$DO_mod_down_sd,
    DO.mod.down.upper = model_output$inst$DO_mod_down_mean + 
      model_output$inst$DO_mod_down_sd
  )
}

# Run function
mod_DO <- return_DO(data_list = data_prepped, model_output = mod_out)

# Plotting
ggplot(data = mod_DO, aes(x = solar.time)) +
  geom_ribbon(aes(ymin = DO.mod.down.lower, ymax = DO.mod.down.upper, 
                  fill = "Modeled"), alpha = 0.4) +
  geom_point(aes(y = DO.obs.down, color = "Observed")) +
  geom_point(aes(y = DO.mod.down, color = "Modeled")) +
  theme_minimal()

# Inspect model fit parameters

# Rhat should be less than 1.2 (to signify convergence)
mod_out$daily$GPP_Rhat
mod_out$daily$ER_Rhat
mod_out$daily$K600_daily_Rhat
mod_out$overall$err_obs_iid_sigma_Rhat
mod_out$overall$lp___Rhat

# Take a look at R2 of O2 fit
mod_out$daily$DO_R2_mean

# Save model results ###########################################################
# Notes: O2 script looks good and runs. Will need adjustments for when modeling more than 1 ply (day) of data. Note that mod_out$inst should probably be trimmed to remove initial lag?
# Next tasks:
# diel n2 model structure
# build in support for multiday datasets
