# Function for extracting relevant daily rhats and r2 from stan model fit object
getdailyfit <- function(modeloutput){
  # Grab stan output from model object
  stanoutput <- modeloutput$stan_out
  
  # Names of Rhat columns in the daily model output
  dailyhat <- names(stanoutput$daily)[str_detect(names(stanoutput$daily), 
                                                 pattern = "Rhat")]
  # Detect R2 columns in daily model output
  dailyR2 <- names(stanoutput$daily)[str_detect(names(stanoutput$daily), 
                                                pattern = "R2_mean")]
  
  # Glue into a single dataframe
  data.frame(date = stanoutput$daily$date,
             stanoutput$daily[dailyhat],
             stanoutput$daily[dailyR2])
}