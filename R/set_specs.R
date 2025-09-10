#' \code{set_specs} Generate a list of model specifications
#' 
#' This function is based on the function of the same name from 
#' `streamMetabolizer`, just simplified to only parameters relavant for 
#' our model type (Bayes, observation error, reaeration depends on modeled gas
#' at downstream station). Default values are the same as `streamMetabolizer`.
#' 
set_specs <- function(
    model_name,
    
    keep_mcmc = TRUE,
    keep_mcmc_data = TRUE,
    
    # Hyperparameters for non-hierarchial GPP, ER
    GPP_daily_mu = 3.1,
    GPP_daily_lower = -Inf,
    GPP_daily_sigma = 6.0,
    #alpha_meanlog = -4.6
    #alpha_sdlog = 0.5,
    #Pmax_mu = 10,
    #Pmax_sigma = 7,
    ER_daily_mu = -7.1,
    ER_daily_upper = Inf,
    ER_daily_sigma = 7.1,
    
    # We're not pooling K600 (we only have 1 day of data)
    # Therefore hyperparameters for non-hierarchial K600
    K600_daily_meanlog = log(12),
    K600_daily_sdlog = 0.05, # not totally sure this is correct
    
    # hyperparameters for error terms
    err_obs_iid_sigma_scale = 0.03,
    
    # Bayes model options
    n_chains = 4,
    n_cores = 4,
    burnin_steps = 500,
    saved_steps = 500,
    thin_steps = 1,
    verbose = FALSE){
  
  # Names of input parameters for model
  params_in <- 
    c("GPP_daily_mu", "GPP_daily_lower", "GPP_daily_sigma",
      "ER_daily_mu", "ER_daily_upper", "ER_daily_sigma",
      "K600_daily_meanlog", "K600_daily_sdlog", "err_obs_iid_sigma_scale")
  
  # Names of output parameters for model
  params_out <- 
    c("GPP", "ER", "DO_R2", "GPP_daily", "ER_daily", "K600_daily",
      "K600_daily_sdlog", "err_obs_iid_sigma")
  
  specsList <- list(
    model_name = model_name,
    
    params_in = params_in,
    params_out = params_out,
    
    keep_mcmc = keep_mcmc,
    keep_mcmc_data = keep_mcmc_data,
    
    GPP_daily_mu = GPP_daily_mu,
    GPP_daily_lower = GPP_daily_lower,
    GPP_daily_sigma = GPP_daily_sigma,
    ER_daily_mu = ER_daily_mu,
    ER_daily_upper = ER_daily_upper,
    ER_daily_sigma = ER_daily_sigma,
    
    K600_daily_meanlog = K600_daily_meanlog,
    K600_daily_sdlog = K600_daily_sdlog, 
    
    err_obs_iid_sigma_scale = err_obs_iid_sigma_scale,
    
    n_chains = n_chains,
    n_cores = n_cores,
    burnin_steps = burnin_steps,
    saved_steps = saved_steps,
    thin_steps = thin_steps,
    verbose = verbose
  )
  return(specsList)
}

