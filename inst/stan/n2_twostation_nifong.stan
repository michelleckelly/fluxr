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
  real K600_daily_meanlog;
  real<lower=0> K600_daily_sdlog;

  // error distributions
  real<lower=0> err_obs_iid_sigma_scale;

  // data dimensions
  int<lower=1> d; // number of dates in time series
  int<lower=1> n24; // max number of observations in 24 hours
  int<lower=1> n; // number of observations per date
  int<lower=1> lag; // lag on date
  real<lower=0> tt; // travel time on date

  // values that change each day
  //vector[d] tt;
  //array[d] int lag;
  //vector[d] n; // number of non-NA observations on each date

  // data values that change throughout day
  vector[d] N2_obs_up[n];
  vector[d] N2_obs_down[n];
  vector[d] N2_sat_up[n];
  vector[d] N2_sat_down[n];
  vector[d] light_mult_N2consume[n];
  vector[d] const_mult_DN[n];
  vector[d] KN2_conv[n];
  vector[d] depth[n];
  vector[d] temp[n];
}

//transformed data{
//  int<lower=1> daylag; // Temp holder for each day's lag value
//  
//  // Grab index start value for day
//  daylag = lag[d] + 1;
//  
//  // Debugging
//  print("daylag = ", daylag);
//}

parameters {
  vector<lower=N2consume_daily_lower>[d] N2consume_daily;
  vector<upper=DN_daily_upper>[d] DN_daily;
  vector<lower=0>[d] K600_daily;

  real<lower=0> err_obs_iid_sigma_scaled;
}

transformed parameters {
  real<lower=0> err_obs_iid_sigma;
  vector[d] N2consume_inst[n];
  vector[d] DN_inst[n];
  vector[d] KN2_inst[n];
  vector[d] N2_mod_down[n];

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
  for(i in 1:n){
    N2_mod_down[i] =
    (
      N2_obs_up[i] +
      (N2consume_inst[i] ./ depth[i]) +
      (DN_inst[i] ./ depth[i] .* tt) +
      KN2_inst[i] .* tt .* (N2_sat_up[i] - N2_obs_up[i] + N2_sat_down[i]) / 2
      ) ./
      ((1 + KN2_inst[i] .* tt / 2));
  }
  
  // Debugging
  //print("N2_mod_down = ", N2_mod_down);
  //print("N2consume_inst = ", N2consume_inst);
  //print("DN_inst = ", DN_inst);
  //print("KN2_inst = ", KN2_inst);
}

model{
  // IID observation error
  for(i in 1:n){
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
  vector[d] err_obs_iid[n];
  vector[d] N2consume;
  vector[d] DN;
  vector[n] N2_obs_down_vec; // temporary, needed to reorganize matrix structure
  vector[n] N2_mod_down_vec; // temporary
  vector[d] N2_R2;
  
  for(i in 1:n){
    err_obs_iid[i] = N2_mod_down[i] - N2_obs_down[i];
  }
  
  // for ply, the below loops 1:d
  for(j in 1:d){
    N2consume[j] = sum(N2consume_inst[1:n, j]) / n24;
    DN[j] = sum(DN_inst[1:n, j]) / n24;
    
    for(i in 1:n){
      N2_mod_down_vec[i] = N2_mod_down[i,j];
      N2_obs_down_vec[i] = N2_obs_down[i,j];
    }
    
    // R2 for N2 = difference between observed and modeled N2 down
  N2_R2[j] = 1 - sum((N2_mod_down_vec - N2_obs_down_vec) .* (N2_mod_down_vec - N2_obs_down_vec)) / sum((N2_obs_down_vec - mean(N2_obs_down_vec)) .* (N2_obs_down_vec - mean(N2_obs_down_vec)));
  }

}
