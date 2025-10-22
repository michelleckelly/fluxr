// This Stan program defines a two-station N2 model,
// as written by Nifong et al. 2020 and available
// at github.com/rlnifong/Denitrification/
// within insst/executables/nn2_model.stan and 
// rewritten following error-handling specifications of 
// b_np_oipi_tr_plrckm.stan from streamMetabolizer
data {
  // parameters of priors on N2 flux
  real N2consume_daily_mu;
  real N2consume_daily_lower;
  real<lower=0> N2consume_daily_sigma;
  real DN_daily_mu;
  real DN_daily_upper;
  real<lower=0> DN_daily_sigma;
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
  row_vector[n] N2_obs_up;
  row_vector[n+lag] N2_obs_down;
  row_vector[n] N2_sat_up;
  row_vector[n+lag] N2_sat_down;
  //row_vector[n] lightfrac;
  row_vector[n] light_mult_N2consume;
  row_vector[n] const_mult_DN;
  row_vector[n] KN2_conv;
  row_vector[n] depth;
  row_vector[n] temp;
}

parameters {
  // NOTE: for ply, these become vectors: vector<lower=0>[d] K600_daily;
  real<lower=N2consume_daily_lower> N2consume_daily;
  real<upper=DN_daily_upper> DN_daily;
  real<lower=0> K600_daily;

  real<lower=0> err_obs_iid_sigma_scaled;
}

transformed parameters {
  real<lower=0> err_obs_iid_sigma;
  row_vector[n] N2consume_inst;
  row_vector[n] DN_inst;
  row_vector[n] KN2_inst;
  row_vector[n+lag] N2_mod_down;

  // Rescale error distribution parameters
  err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled;

  // Model N2 time series
  // * two-station
  // * observation error
  // * reaeration depends on N2_mod_down
  // * N2consume is linked to light

  // Calculate individual process rates
  for (i in 1:n) {
    N2consume_inst[i] = N2consume_daily .* light_mult_N2consume[i];
    DN_inst[i] = DN_daily .* const_mult_DN[i];
    KN2_inst[i] = K600_daily .* KN2_conv[i];
  }
  
  // N2 model
  N2_mod_down[1:lag] = N2_obs_down[1:lag];
  for(i in 1:n){
    N2_mod_down[i+lag] =
    (
      N2_obs_up[i] +
      (N2consume_inst[i] ./ depth[i]) +
      (DN_inst[i] ./ depth[i] * tt) +
      KN2_inst[i] * tt .* (N2_sat_up[i] - N2_obs_up[i] + N2_sat_down[i+lag]) / 2
      ) ./
      ((1 + KN2_inst[i] * tt / 2));
  }
  
  // Debugging
  //print("N2_mod_down = ", N2_mod_down);
  //print("N2consume_inst = ", N2consume_inst);
  //print("DN_inst = ", DN_inst);
  //print("KN2_inst = ", KN2_inst);
}

model{
  // IID observation error
  for(i in 1:(n+lag)){
    N2_obs_down[i] ~ normal(N2_mod_down[i], err_obs_iid_sigma);
  }
  // SD of observation error
  err_obs_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  N2consume_daily ~ normal(N2consume_daily_mu, N2consume_daily_sigma);
  DN_daily ~ normal(DN_daily_mu, DN_daily_sigma);
  K600_daily ~ lognormal(K600_daily_meanlog, K600_daily_sdlog);
}

generated quantities{
  row_vector[n+lag] err_obs_iid;
  real N2consume;
  real DN;
  real N2_R2;

  for(i in 1:(n+lag)){
    err_obs_iid[i] = N2_mod_down[i] - N2_obs_down[i];
  }
  
  N2consume = sum(N2consume_inst[1:n]) / n;
  DN = sum(DN_inst[1:n]) / n;
  
  // R2 = difference between observed and modeled N2 down
  N2_R2 = 1 - sum((N2_mod_down - N2_obs_down) .* (N2_mod_down - N2_obs_down)) / sum((N2_obs_down - mean(N2_obs_down)) .* (N2_obs_down - mean(N2_obs_down)));

}

