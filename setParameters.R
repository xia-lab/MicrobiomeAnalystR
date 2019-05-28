setParameters <- function(file_compressed = TRUE, # logic TURE or FALSE, compressed format gz, bz2,
                          OS_is_windows = FALSE, # operation system is windows, linux or mac;
                          ... ){
  parameters <- list(file_compressed = file_compressed,
                     OS_is_windows = OS_is_windows);
  return(parameters);
}
