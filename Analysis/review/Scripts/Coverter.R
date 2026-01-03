convert_edf <- function(data_dir = "data") {
  library(processx)
  # Get all subdirectories
  all_dirs <- list.dirs(data_dir, recursive = TRUE, full.names = TRUE)
  
  for (dir in all_dirs) {
    edf_files <- list.files(dir, pattern = "\\.EDF$", ignore.case = TRUE, full.names = TRUE)
    asc_files <- list.files(dir, pattern = "\\.asc$", ignore.case = TRUE, full.names = TRUE)
    
    if (length(edf_files) > 0 && length(asc_files) == 0) {
      message("Converting in: ", dir)
      for (edf in edf_files) {
        # Use system call to convert
        cmd <- sprintf('"%s" -y "%s"', "edf2asc", edf)
        system(cmd)
      }
    } else if (length(edf_files) > 0 && length(asc_files) > 0) {
      message("Skipping (ASC already exists): ", dir)
    }
  }
}






