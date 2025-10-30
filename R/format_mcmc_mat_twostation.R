#' \code{format_mcmc_mat_twostation} Format MCMC output into a one-row data.frame
#'
#' @param mcmc_mat matrix as extracted from Stan
#' @param names_params character vector of the names of the parameters
#' @param names_stats character vector of the names of the statistics
#' @import dplyr
#' @keywords internal
format_mcmc_mat_twostation <- function(mcmc_mat, #data_list_d, data_list_n,
                                       keep_mcmc, runmcmc_out) {
  
  # This function is basically entirely ripped from 
  # streamMetabolizer::format_mcmc_mat_nosplit, just with some changes for the 
  # N2 terms. streamMetabolizer's index reformatting steps are very clever
  par_homes <-
    list(
      overall = c(
        "err_obs_iid_sigma", "err_obs_iid_sigma_scaled", "lp__"
      ),
      daily = c(
        "K600_daily", 
        "GPP", "ER", "GPP_daily", "ER_daily", "DO_R2",
        "N2consume", "DN", "N2consume_daily", "DN_daily", "N2_R2"
        ),
      inst = c(
        "err_obs_iid",
        "DO_mod_down", "GPP_inst", "ER_inst", "KO2_inst", 
        "N2_mod_down", "N2consume_inst", "DN_inst", "KN2_inst")
    )
  
  # Declare dplyr variables
  stat <- val <- . <- rowname <- variable <- varstat <-
    indexstr <- date_index <- time_index <- index <- '.dplyr_var'
  
  # Using rownames of mcmc_mat, determine which data.frames to create and 
  # which params should be included in each
  var_table <- table(gsub("\\[[[:digit:]|,]+\\]", "", rownames(mcmc_mat)))
  par_dfs <- sapply(names(var_table), function(parname){
    home <- names(par_homes)[which(sapply(par_homes, function(vd) parname %in% vd))]
    if(length(home) == 0){
      home <- var_table[[parname]]
    }
    return(home)
  })
  all_dims <- lapply(setNames(nm = unique(par_dfs)), 
                     function(upd) {names(par_dfs)[which(par_dfs == upd)]}
                     )
  
  # Replace percentage symbol
  colnames(mcmc_mat) <- gsub("%", "pct", colnames(mcmc_mat))
  # Create data.frame with vars in columns and indicies in rows
  mcmc_out <- lapply(setNames(nm = names(all_dims)), function(dfname){
    df_params <- all_dims[[dfname]]
    dim_rows <- sort(do.call(c, lapply(df_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat)))))
    row_order <- names(sort(sapply(df_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat))[1])))
    varstat_order <- paste0(rep(row_order, each = ncol(mcmc_mat)), "_", rep(colnames(mcmc_mat), times = length(row_order)))
    par_dims <- sapply(df_params, function(dp) length(grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat))))
    
    tibble::as_tibble(mcmc_mat[dim_rows, , drop = FALSE]) %>%
      mutate(rowname = rownames(mcmc_mat[dim_rows, , drop = FALSE])) %>%
      select(rowname, everything()) %>%
      gather(stat, value = val, 2:ncol(.)) %>%
      mutate(variable = gsub("\\[[[:digit:]|,]+\\]", "", rowname),
             indexstr = if(1 %in% par_dims) "1" else sapply(strsplit(rowname, "\\[|\\]"), "[[", 2),
             # parse the index number
             index = 
               if(dfname %in% c("daily", "inst") || any(grepl(",", indexstr))) {
                 indexstr
               } else {
                 as.numeric(indexstr)
               },
             date_index = 
               if(dfname == "daily"){
                 as.numeric(indexstr)
               } else if(dfname == "inst" || any(grepl(",", indexstr))) {
                 sapply(strsplit(indexstr, ","), function(ind) as.numeric(ind[2]))
               } else{
                 NA
               },
             time_index = 
               if(dfname == "inst" || any(grepl(",", indexstr))) {
                 sapply(strsplit(indexstr, ","), function(ind) as.numeric(ind[1]))
               } else{
                 NA
               },
             # this is to order the statistics for each variable
             varstat = ordered(paste0(variable, "_", stat), varstat_order)) %>%
      select(date_index, time_index, index, varstat, val) %>%
      arrange(date_index, time_index, index) %>%
      spread(varstat, val)
  })
  
  # add in the model object if user requested
  if(keep_mcmc == TRUE){
    mcmc_out <- c(mcmc_out, list(mcmc_fit = runmcmc_out))
  }
  
  return(mcmc_out)
}
