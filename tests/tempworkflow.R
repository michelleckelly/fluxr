# Pull together an example multi-day dataset
library(tidyverse)
# Load in raw metabolism parameters

# Load in clean metabolism parameters
martha <- read_rds("MetabolismParameters_Clean_MART_2019-01_2019-12.rds")

# Save as csv file
filt <- 
  martha$cleanData %>%
  filter(solarTime >= ymd_hms("2019-05-01 00:00:00") &
           solarTime <= ymd_hms("2019-10-31 23:59:59"))

write_csv(filt, "multidayO2.csv")

# editing data files
#multidayO2 <- read_csv("inst/extdata/multidayO2.csv")
#usethis::use_data(multidayO2, overwrite = TRUE)

### Workflow
data("multidayO2")

# Filter for testing
data_all <- 
  multidayO2 %>%
  filter(solar.time <= ymd_hms("2019-05-11 00:00:00"))

# Modeling settings
modname <- "o2_twostation"
mod_specs <- 
  set_specs(model_name = modname, verbose = TRUE,
            K600_daily_meanlog = 1, K600_daily_sdlog = 0.05,
            n_cores = 8, n_chains = 8)

upname <- "S1"
downname <- "S2"

## new wrapper function for prepping data format and executing model

run_twostation <- function(data_all, mod_specs, upname, downname){
  # Prep format of data
  data_prepped <- 
    prepdata_bayes_twostation(data = data_all, specs = mod_specs,
                              up.name = upname, down.name = downname)
  
  # execute the bayes function inside an error handler
  # Possibly use error handler to capture the undefined values warnings?
  bayes_allday <- runstan_bayes_twostation(data_list = data_prepped, 
                                           specs = mod_specs)
  
  # NOTE: STOPPED HERE
  # Match date and time info to indicies
  # see L370 of metab_bayes.R
  
  
#}