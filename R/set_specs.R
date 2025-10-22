#' \code{set_specs} Generate a list of model specifications
#' 
#' This function is based on the function of the same name from 
#' `streamMetabolizer`, just simplified to only parameters relavant for 
#' our model type (Bayes, observation error, reaeration depends on modeled gas
#' at downstream station).  
#' 
#' Default values for GPP, ER, and K600 are the same as `streamMetabolizer`. 
#' Default values for N2consume and DN are the same as those in Nifong et al. 
#' 2020.
#' 
set_specs <- function(
    model_name,
    
    keep_mcmc = TRUE,
    keep_mcmc_data = TRUE,
    
    # Default hyperparameters for non-hierarchial GPP, ER
    GPP_daily_mu = 3.1,
    GPP_daily_lower = -Inf,
    GPP_daily_sigma = 6.0,
    ER_daily_mu = -7.1,
    ER_daily_upper = Inf,
    ER_daily_sigma = 7.1,
    
    # Default hyperparameters for N2consume, DN
    # As used in Nifong et al 2020
    N2consume_daily_mu = -1,
    N2consume_daily_lower = -Inf,
    N2consume_daily_sigma = 5,
    DN_daily_mu = 1,
    DN_daily_upper = Inf,
    DN_daily_sigma = 5,
    
    # We're not pooling K600 (we only have 1 day of data)
    # Therefore hyperparameters for non-hierarchial K600
    K600_daily_meanlog = log(12),
    K600_daily_sdlog = 0.05,
    
    # hyperparameters for error terms
    err_obs_iid_sigma_scale = 0.03,
    
    # Bayes model options
    n_chains = 4,
    n_cores = 4,
    burnin_steps = 500,
    saved_steps = 500,
    thin_steps = 1,
    verbose = FALSE){
  
  # Set common input parameters for all models
  params_in <- 
    c("K600_daily_meanlog", "K600_daily_sdlog", "err_obs_iid_sigma_scale")
  
  # Set common output parameters for all models
  params_out <-
    c("K600_daily", "err_obs_iid_sigma", "err_obs_iid_sigma_scaled",
      "err_obs_iid")
  
  # Set common contents for specification list
  specsList <- list(
    model_name = model_name,
    
    keep_mcmc = keep_mcmc,
    keep_mcmc_data = keep_mcmc_data,
    
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
  
  # Set additional parameters for O2 models & additional outputs for specsList
  if(model_name == "o2_twostation"){
    
    params_in <- 
      c(params_in, 
        "GPP_daily_mu", "GPP_daily_lower", "GPP_daily_sigma",
        "ER_daily_mu", "ER_daily_upper", "ER_daily_sigma")
    
    params_out <-
      c(params_out,
        "GPP", "ER", "DO_R2", "GPP_daily", "ER_daily", 
        "KO2_inst", "DO_mod_down", "GPP_inst", "ER_inst")
    
    # Add model-specific parameters to specification list
    specsList$params_in <- params_in
    specsList$params_out <- params_out
    
    specsList$GPP_daily_mu <- GPP_daily_mu
    specsList$GPP_daily_lower <- GPP_daily_lower
    specsList$GPP_daily_sigma <- GPP_daily_sigma
    specsList$ER_daily_mu <- ER_daily_mu
    specsList$ER_daily_upper <- ER_daily_upper
    specsList$ER_daily_sigma <- ER_daily_sigma
  }
  
  # Set additional parameters for N2 models, output specsList
  if(model_name == "n2_twostation_nifong"){
    
    params_in <- 
      c(params_in, 
        "N2consume_daily_mu", "N2consume_daily_lower", "N2consume_daily_sigma",
        "DN_daily_mu", "DN_daily_upper", "DN_daily_sigma")
    
    params_out <-
      c(params_out,
        "N2consume", "DN", "N2_R2", "N2consume_daily", "DN_daily", 
        "KN2_inst", "N2_mod_down", "N2consume_inst", "DN_inst")
    
    specsList$params_in <- params_in
    specsList$params_out <- params_out
    
    specsList$N2consume_daily_mu <- N2consume_daily_mu
    specsList$N2consume_daily_lower <- N2consume_daily_lower
    specsList$N2consume_daily_sigma <- N2consume_daily_sigma
    specsList$DN_daily_mu <- DN_daily_mu
    specsList$DN_daily_upper <- DN_daily_upper
    specsList$DN_daily_sigma <- DN_daily_sigma
  }
  
  return(specsList)
}

