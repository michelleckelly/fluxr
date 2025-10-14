// This Stan program defines a two-station O2 model,
// as written by Nifong et al. ... and available
// at github.com/rlnifong/Denitrification/
// within insst/executables/nn2_model.stan
//
// Written following specifications of b_np_oipi_tr_plrckm.stan from
// streamMetabolizer
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
  //int<lower=1> d; // number of dates in time series
  int<lower=1> n24; // number of observations in first 24 hours, per date
  int<lower=1> n; // number of observations per date

  // data values that are fixed each day
  real<lower=0> tt;
  int<lower=0> lag;

  // data values that change throughout day
  // NOTE When you add plys, you'll need to modify these to be vectors rather 
  // than row_vector, i.e. vector[d] DO_obs_up[n];
  row_vector[n] DO_obs_up;
  row_vector[n+lag] DO_obs_down;
  row_vector[n] DO_sat_up;
  row_vector[n+lag] DO_sat_down;
  //row_vector[n] lightfrac;
  row_vector[n] light_mult_GPP;
  row_vector[n] const_mult_ER;
  row_vector[n] KO2_conv;
  row_vector[n] depth;
  row_vector[n] temp;
}

parameters {
  // NOTE: for ply, these become vectors: vector<lower=0>[d] K600_daily;
  real<lower=GPP_daily_lower> GPP_daily;
  real<upper=ER_daily_upper> ER_daily;
  real<lower=0> K600_daily;

  real<lower=0> err_obs_iid_sigma_scaled;
  //real<lower=0> err_proc_iid_sigma_scaled;
  
  //row_vector[n+lag] DO_mod_down;
}

transformed parameters {
  real<lower=0> err_obs_iid_sigma;
  //row_vector[n+lag] DO_mod_down_partial_sigma;
  row_vector[n] GPP_inst;
  row_vector[n] ER_inst;
  row_vector[n] KO2_inst;
  row_vector[n+lag] DO_mod_down;

  // Rescale error distribution parameters
  err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled;

  // Model DO time series
  // * two-station
  // * observation error
  // (two-station is already a proccess error model because we're working with
  // the difference between two known quantities of O2)
  // * reaeration depends on DO_mod_down

  // Calculate individual process rates
  for (i in 1:n) {
    GPP_inst[i] = GPP_daily .* light_mult_GPP[i];
    ER_inst[i] = ER_daily .* const_mult_ER[i];
    KO2_inst[i] = K600_daily .* KO2_conv[i];
  }
  
  // DO model
  // Set DO_mod_down values prior to i+lag = DO_obs_down[i:lag] to avoid
  // NA values in DO_mod_down[i:lag]
  DO_mod_down[1:lag] = DO_obs_down[1:lag];
  for(i in 1:n){
    DO_mod_down[i+lag] =
    (
      DO_obs_up[i] +
      (GPP_inst[i] ./ depth[i]) +
      (ER_inst[i] ./ depth[i] * tt) +
      KO2_inst[i] * tt .* (DO_sat_up[i] - DO_obs_up[i] + DO_sat_down[i+lag]) / 2
      ) ./
      ((1 + KO2_inst[i] * tt / 2));
  }
  
  // Debugging
  //print("DO_mod_down = ", DO_mod_down);
  //print("GPP_inst = ", GPP_inst);
  //print("ER_inst = ", ER_inst);
  //print("KO2_inst = ", KO2_inst);
}
model{
  // IID observation error
  for(i in 1:(n+lag)){
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
  // ply: below become vectors vector[d] err_obs_iid[n]; vector[d] GPP;
  row_vector[n+lag] err_obs_iid;
  real GPP;
  real ER;
  real DO_R2;

  for(i in 1:(n+lag)){
    err_obs_iid[i] = DO_mod_down[i] - DO_obs_down[i];
  }
  // for ply, the below loops 1:d
  // these were originally sum(GPP_inst[1:n24]) / n24; because of assumption
  // that we've checked for full days. divide by n for now
  GPP = sum(GPP_inst[1:n]) / n;
  ER = sum(ER_inst[1:n]) / n;
  
  // R2 for DO = difference between observed and modeled DO down
  DO_R2 = 1 - sum((DO_mod_down - DO_obs_down) .* (DO_mod_down - DO_obs_down)) / sum((DO_obs_down - mean(DO_obs_down)) .* (DO_obs_down - mean(DO_obs_down)));

 // DO_mod_down_vec = DO_mod_down;
}
