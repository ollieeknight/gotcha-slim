.onAttach <- function(libname, pkgname) {
  library(future, quietly = TRUE)
  library(future.apply, quietly = TRUE)
  library(ShortRead, quietly = TRUE)
  library(fs, quietly = TRUE)
  library(cli, quietly = TRUE)
  library(stringr, quietly = TRUE)
  library(progressr, quietly = TRUE)
  library(data.table, quietly = TRUE)
  library(Biostrings, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(data.table, quietly = TRUE)


  packageStartupMessage(
    "
░██████╗░░█████╗░████████╗░█████╗░██╗░░██╗░█████╗░  ░░░░░░  ░██████╗██╗░░░░░██╗███╗░░░███╗
██╔════╝░██╔══██╗╚══██╔══╝██╔══██╗██║░░██║██╔══██╗  ░░░░░░  ██╔════╝██║░░░░░██║████╗░████║
██║░░██╗░██║░░██║░░░██║░░░██║░░╚═╝███████║███████║  █████╗  ╚█████╗░██║░░░░░██║██╔████╔██║
██║░░╚██╗██║░░██║░░░██║░░░██║░░██╗██╔══██║██╔══██║  ╚════╝  ░╚═══██╗██║░░░░░██║██║╚██╔╝██║
╚██████╔╝╚█████╔╝░░░██║░░░╚█████╔╝██║░░██║██║░░██║  ░░░░░░  ██████╔╝███████╗██║██║░╚═╝░██║
░╚═════╝░░╚════╝░░░░╚═╝░░░░╚════╝░╚═╝░░╚═╝╚═╝░░╚═╝  ░░░░░░  ╚═════╝░╚══════╝╚═╝╚═╝░░░░░╚═╝
"
  )
  
  # Print the warning message
  cli::cli_alert_warning("To check the Python environment and required packages, run check_gotcha_slim_environment() after setting r-reticulate as your reticulate Python environment.")
}
