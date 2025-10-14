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
            K600_daily_meanlog = 1, K600_daily_sdlog = 0.05)

# Pass data and specs to prepping function #####################################
data_prepped <- 
  prepdata_bayes_twostation(data = data_all, ply_date = "2024-09-19", 
                            specs = mod_specs)

# FOR TESTING
data_list <- data_prepped

# Run model ####################################################################
#runstan_bayes(data_list = data_prepped, specs = mod_specs)

#runstan_bayes <- function(data_list, specs, ){
  # Pull out model name
  model_name <- data_list$model_name
  
  # Pull settings from specs list
  ncores <- data_list$n_cores
  n_chains <- data_list$n_chains
  verbose <- data_list$verbose
  warmup <- data_list$burnin_steps
  iter <- data_list$saved_steps + warmup
  thin_steps <- data_list$thin_steps
  params_out <- data_list$params_out
  keep_mcmc <- data_list$keep_mcmc
  
  # Determine number of cores to use
  tot_cores <- parallel::detectCores()
  # Compare to number of cores specified in set_specs()
  n_cores <- min(tot_cores, ncores)
  
  # Tell user how many cores Stan is going to use so they can adjust settings
  # if they want
  if(verbose){message(paste0("\nStan requesting ", n_chains, " chains on ", 
                             n_cores, " of ", tot_cores, " available cores."))}
  
  # Get filepath to compiled & uncompiled stan model
  model_path <- paste0("inst/stan/", model_name, ".stan")
  mobj_path <- paste0("inst/stan/", model_name, ".rds")
  
  # Check for existence of compiled model. If model doesn't exist, compile it.
  if(file.exists(mobj_path)){
    if(verbose){message(paste0("\nUsing previously compiled ", model_name, 
                               " model: ", mobj_path))}
  } else{
    # Start clock for compilation time
    start.time <- Sys.time()
    # If file doesn't exist, compile model
    if(verbose){message(paste0("\nCompiling ", model_name, ".stan model..."))}
    stan_model(file = model_path, model_name = model_name, auto_write = TRUE,
               warn_pedantic = TRUE, warn_uninitialized = TRUE)
    # Stop clock
    end.time <- Sys.time()
    
    # Verify that compilation was successful & present done message
    if(file.exists(mobj_path)){
      elapsed <- round(as.numeric(end.time - start.time, 
                                  units = "hours") * 60 * 60,
                       digits = 2)
      if(verbose){message(paste0("\nDone! Compile time was ", elapsed, 
                                 " seconds."))}
    } else{
      message("ERROR: Stan model was not compiled successfully. Try restarting your R session.")
    }
  }
  
  # streamMetabolizer has the following garbage collection line with the 
  # comment "this humble line saves us from many horrible R crashes" which
  # sounds good so replicating that here
  gc()
  
  # Read compiled Stan file into memory
  stan_mobj <- readRDS(mobj_path)
  
  # Run Stan
  if(verbose){message("\nSampling Stan model")}
  # Note whatever log files are in the temp directory prior to model run
  oldlogfiles <- 
    normalizePath(file.path(tempdir(), grep("_StanProgress.txt", 
                                            dir(tempdir()), value = TRUE)))
  
  start.time <- Sys.time()
  
  consolelog <- capture.output(
    runstan_out <- rstan::sampling(
      object = stan_mobj,
      data = data_list, 
      pars = params_out, 
      include = TRUE,
      chains = n_chains,
      warmup = warmup,
      iter = iter, 
      thin = thin_steps,
      init = "random",
      verbose = verbose,
      open_progress = FALSE,
      cores = n_cores),
    split = verbose)
  
  # streamMetabolizer puts a breakpoint here when debugging
  
  # View output
  show(runstan_out)
  # Plot estimates
  rstan::plot(runstan_out)
  # Display traceplot of chains
  traceplot(runstan_out)
  # Density plot
  stan_dens(runstan_out)
  
  rstan::extract(runstan_out, "DO_mod_down")
  runstan_out@model_pars
  
  # Extract list of posterior distribution draws
  draws_list <- rstan::extract(runstan_out)
  as.matrix(runstan_out)
  summary(runstan_out)
  
  # Check out what's stored in runstan_out
  names(runstan_out)
  runstan_out@model_name # model name
  runstan_out@model_pars # list of modeled parameters
  runstan_out@inits # model outputs for each chain
  runstan_out@date # modeling date
  
  
  end.time <- Sys.time()
  elapsed <- round(as.numeric(end.time - start.time, 
                              units = "hours") * 60, digits = 2)
  
  # Check for failed model run
  if(runstan_out@mode == 2L){
    # Keep mcmc for failed model run
    stan_out <- NULL
    warning(capture.output(print(runstan_out)))
  } else {
    # Assuming here that we're just running 1-day models
    stan_mat <- rstan::summary(runstan_out)$summary
    # Pull modeled parameter (row) names
    names_params <- rep(gsub("\\[1\\]", "", rownames(stan_mat)), 
                        each = ncol(stan_mat))
    # Pull stats (column) names
    names_stats <- rep(gsub("%", "pct", colnames(stan_mat)), 
                       each = nrow(stan_mat))
    # Format stan output
    stan_out <- format_mcmc_mat(mcmc_mat = stan_mat,# names_params, names_stats, 
                                keep_mcmc, runmcmc_out = runstan_out,
                                data_list)
    
  }
  
  # Attach contents of most recent logfile (which will be for this model)
  newlogfiles <- 
    normalizePath(file.path(tempdir(), grep("_StanProgress.txt", 
                                            dir(tempdir()), value = TRUE)))
  logfile <- setdiff(newlogfiles, oldlogfiles)
  log <- 
    if(length(logfile) > 0){
      readLines(logfile)
      } else{consolelog}
  stan_out <- c(stan_out, c(
    list(log = log),
    if(exists("compile_log")){list(compile_log = compile_log)},
    list(compile_time_min = elapsed)))
  
  return(stan_out)
  
  # Double check that stan_out will look like the streammetabolizer model output (see vignette)
  
  # Q: how to make sure modeled O2 is output also?
#}
