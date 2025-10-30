// This Stan program defines a two-station O2 model,
// as written by Nifong et al. 2020 and available
// at github.com/rlnifong/Denitrification/
// within insst/executables/nn2_model.stan and 
// rewritten following error-handling specifications of 
// b_np_oipi_tr_plrckm.stan from streamMetabolizer
data {
  // parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_lower;
  real<lower=0> GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_upper;
  real<lower=0> ER_daily_sigma;
  real<lower=0> K600_daily_meanlog;
  real<lower=0> K600_daily_sdlog;

  // error distributions
  real<lower=0> err_obs_iid_sigma_scale;

  // data dimensions
  int<lower=1> d; // number of dates in time series
  int<lower=1> n24; // max number of observations in 24 hours
  //int<lower=1> n; // number of observations per date
  
  // values that change each day
  vector[d] tt;
  array[d] int lag;
  vector[d] n; // number of non-NA observations on each date

  // data values that change throughout day
  vector[d] DO_obs_up[n24];
  vector[d] DO_obs_down[n24];
  vector[d] DO_sat_up[n24];
  vector[d] DO_sat_down[n24];
  vector[d] light_mult_GPP[n24];
  vector[d] const_mult_ER[n24];
  vector[d] KO2_conv[n24];
  vector[d] depth[n24];
  vector[d] temp[n24];
  
}
transformed data{
  int<lower=1> daylag; // Temp holder for each day's lag value
  
  // Grab index start value for day
  daylag = lag[d] + 1;
  
  // Debugging
  print("daylag = ", daylag);
}

parameters {
  vector<lower=GPP_daily_lower>[d] GPP_daily;
  vector<upper=ER_daily_upper>[d] ER_daily;
  vector<lower=0>[d] K600_daily;

  real<lower=0> err_obs_iid_sigma_scaled;
}

transformed parameters {
  real<lower=0> err_obs_iid_sigma;
  vector[d] GPP_inst[n24];
  vector[d] ER_inst[n24];
  vector[d] KO2_inst[n24];
  vector[d] DO_mod_down[n24];
  
  // Rescale error distribution parameters
  err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled;
  
  // Model DO time series
  // * two-station
  // * observation error
  // (two-station is already a proccess error model because we're working with
  // the difference between two known quantities of O2)
  // * reaeration depends on DO_mod_down

  // Calculate individual process rates
  for (i in daylag:n24) {
    GPP_inst[i] = GPP_daily .* light_mult_GPP[i];
    ER_inst[i] = ER_daily .* const_mult_ER[i];
    KO2_inst[i] = K600_daily .* KO2_conv[i];
  }
  
  // DO model
  
  // Set DO_mod_down values prior to i+lag = DO_obs_down[i:lag] to avoid
  // NA values in DO_mod_down[i:lag]
  //DO_mod_down[1+lag] = DO_obs_down[1+lag];
  for(i in daylag:n24){
    DO_mod_down[i] =
    (
      DO_obs_up[i] +
      (GPP_inst[i] ./ depth[i]) +
      (ER_inst[i] ./ depth[i] .* tt) +
      KO2_inst[i] .* tt .* (DO_sat_up[i] - DO_obs_up[i] + DO_sat_down[i]) / 2
      ) ./
      ((1 + KO2_inst[i] .* tt / 2));
  }
  
  // Debugging
  //print("daylag = ", daylag);
  //print("DO_mod_down = ", DO_mod_down);
  //print("GPP_inst = ", GPP_inst);
  //print("ER_inst = ", ER_inst);
  //print("KO2_inst = ", KO2_inst);
}
model{
  // IID observation error
  for(i in daylag:n24){
    DO_obs_down[i] ~ normal(DO_mod_down[i], err_obs_iid_sigma);
  }
  // SD of observation error
  err_obs_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ lognormal(K600_daily_meanlog, K600_daily_sdlog);
}
generated quantities{
  vector[d] err_obs_iid[n24];
  vector[d] GPP;
  vector[d] ER;
  vector[n24] DO_obs_down_vec; // temporary, needed to reorganize matrix structure
  vector[n24] DO_mod_down_vec; // temporary
  vector[d] DO_R2;

  for(i in daylag:n24){
    err_obs_iid[i] = DO_mod_down[i] - DO_obs_down[i];
  }
  // for ply, the below loops 1:d
  // these were originally sum(GPP_inst[1:n24]) / n24; because of assumption
  // that we've checked for full days. divide by n for now
  for(j in 1:d){
    GPP[j] = sum(GPP_inst[1:n24, j]) / n24;
    ER[j] = sum(ER_inst[1:n24, j]) / n24;
    
    for(i in daylag:n24){
      DO_mod_down_vec[i] = DO_mod_down[i,j];
      DO_obs_down_vec[i] = DO_obs_down[i,j];
    }
    
    // R2 for DO = difference between observed and modeled DO down
  DO_R2[j] = 1 - sum((DO_mod_down_vec - DO_obs_down_vec) .* (DO_mod_down_vec - DO_obs_down_vec)) / sum((DO_obs_down_vec - mean(DO_obs_down_vec)) .* (DO_obs_down_vec - mean(DO_obs_down_vec)));
  }
  
}
