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
  rename(light = light.calc)

# Add date column
data_all$date <- date(data_all$solar.time)

# Modeling settings ############################################################
modname <- "n2_twostation_nifong" # N2
#modname <- "o2_twostation" # O2

mod_specs <- 
  set_specs(model_name = modname, verbose = TRUE,
            K600_daily_meanlog = 1, K600_daily_sdlog = 0.05,
            n_cores = 8, n_chains = 8)

# Prep data ####################################################################
data_prepped <- 
  prepdata_bayes_twostation(data = data_all, ply_date = "2024-09-19", 
                            specs = mod_specs)

# Run model ####################################################################
mod_out <- 
  runstan_bayes_twostation(data_list = data_prepped, 
                           specs = mod_specs)

# Evaluate fit and predictions #################################################
# View daily fit statistics 
getdailyfit(mod_out)

# View predicted rates
# Flux units are g-[O2 N2] m-2 d-1
mod_out$stan_out$daily %>%
  select(date, matches("daily_mean|daily_sd"))

# For mg-N m-2 h-1: value / 28 * 2 * 14 * 1000 / 24

mod_out$stan_out$daily$DN_daily_mean * 1000 / 24
mod_out$stan_out$daily$N2consume_daily_mean * 1000 / 24


# Plot modeled vs observed gas
if(modname == "n2_twostation_nifong"){
  ggplot(data = mod_out$stan_out$conc, 
         aes(x = solar.time)) +
    geom_ribbon(aes(ymin = N2.mod.down.lower, ymax = N2.mod.down.upper, 
                    fill = "Modeled"), alpha = 0.4) +
    geom_point(aes(y = N2.obs.down, color = "Observed")) +
    geom_point(aes(y = N2.mod.down, color = "Modeled")) +
    theme_minimal()
  } else{
    ggplot(data = mod_out$stan_out$conc, 
           aes(x = solar.time)) +
      geom_ribbon(aes(ymin = DO.mod.down.lower, ymax = DO.mod.down.upper, 
                      fill = "Modeled"), alpha = 0.4) +
      geom_point(aes(y = DO.obs.down, color = "Observed")) +
      geom_point(aes(y = DO.mod.down, color = "Modeled")) +
      theme_minimal()
  }

# Save model results ###########################################################
#filename <- ""

# Next tasks:
# build in support for multiday datasets
