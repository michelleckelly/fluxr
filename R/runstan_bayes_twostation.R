#' \code{runstan_bayes_twostation}  Fit two-station Bayesian metabolism model
#'
#' @param data_list list of inputs to Stan model, as returned by `prepdata_bayes_twostation()`
#' @param specs list of model specifications, as returned by `set_specs()`
#' @import dplyr
runstan_bayes_twostation <- function(data_list, specs){
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
    rstan::stan_model(file = model_path, model_name = model_name, 
                      auto_write = TRUE, warn_pedantic = TRUE, 
                      warn_uninitialized = TRUE)
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
  
  # Detect if run was successful
  end.time <- Sys.time()
  elapsed <- end.time - start.time
  
  # streamMetabolizer puts a breakpoint here when debugging 
  
  # Check for a failed model run
  if(runstan_out@mode == 2L){
    # Keep mcmc for failed model run
    stan_out <- NULL
    warning(capture.output(print(runstan_out)))
  } else {
    stan_mat <- rstan::summary(runstan_out)$summary
    # Format stan output
    stan_out <- format_mcmc_mat_twostation(stan_mat, keep_mcmc, runstan_out)
    
    ## NOTE -- ADD HERE ##
    # Add dataframe of predicted and actual gas concentrations
    #stan_out$conc <- 
     # return_gas(modname = model_name, 
      #           instresults = stan_out$inst, 
       #          inputdata = data_list)
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
    list(compile_time = elapsed)))
  
  return(stan_out)
}
