#' Look for a model file, as written by streamMetabolizer
#' 
#' Looks first in the models folder of the fluxr package, second 
#' along the relative or absolute file path given by model_name
#' 
#' @param model_name a model file in the 'models' folder of the 
#'   fluxr package or a relative or absolute file path of a model 
#'   file
#' @return a file path if the file exists or an error otherwise
#' @keywords internal
locate_filename_fluxr <- function(model_name) {
  model_name <- paste0(model_name, ".rds")
  
  package_dir <- system.file("stan", package="fluxr")
  package_path <- system.file(paste0("stan/", model_name), package="fluxr")
  other_path <- model_name
  
  if(file.exists(package_path)) return(package_path)
  if(file.exists(other_path)) return(other_path)
  
  #message("could not locate the model file at ", file.path(package_dir, model_name), " or ", other_path) 
  return(paste0(package_dir, "/", other_path))
}