setup_explore <- function(options){
  imfil_explore <- options$explore
  explore_function <- c()
  explore_data_flag <- c()
  explore_data <- c()
  if (imfil_explore == 1){
    explore_function   <- options$explore_function
    explore_data_flag  <- options$explore_data_flag
    if (explore_data_flag == 1) explore_data <- options$explore_data
  }
  list(explore_function = explore_function, explore_data_flag = explore_data_flag, 
       explore_data = explore_data)
}