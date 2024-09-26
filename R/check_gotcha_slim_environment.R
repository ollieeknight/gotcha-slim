check_gotcha_slim_environment <- function(conda_env_name = "r-reticulate") {
  
  cli::cli_h1("Checking reticulate environment")
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli::cli_abort("The 'reticulate' package is required but not installed")
  }

  # Get the Conda environment path
  conda_envs <- reticulate::conda_list()
  conda_env_path <- conda_envs$python[conda_envs$name == conda_env_name]

  if (length(conda_env_path) == 0) {
    cli::cli_abort(
      "{.red Conda environment not found.} Please create the Conda environment '{.envvar {conda_env_name}}' with the required packages."
    )
  }

  reticulate::use_condaenv(conda_env_name, required = TRUE)

  required_packages <-  c("pandas", "matplotlib", "seaborn", "numpy", "sklearn")

  installed_packages <- tryCatch(
    {
      sapply(required_packages, function(pkg) {
        reticulate::py_module_available(pkg)
      })
    },
    error = function(e) {
      cli::cli_alert_danger("Error checking installed Python packages")
      NULL
    }
  )

  missing_packages <- required_packages[!installed_packages]

  if (length(missing_packages) > 0) {
    cli::cli_abort(
      "{.red The following required Python packages are not installed:} {.pkg {paste(missing_packages, collapse = ', ')}}. \n",
      "Hint: Install the missing packages using `conda install {paste(missing_packages, collapse = ' ')} --name {conda_env_name}`."
    )
  } else {
    cli::cli_alert_success("All required Python packages are already installed.")
  }

  writeLines(conda_env_path, con = paste0(find.package('gotchaslim'), '/last_used_environment.txt'))

  # Download the file
  download.file(
    "https://raw.githubusercontent.com/ollieeknight/gotcha-slim/main/inst/python/gotcha_labeling.py",
    paste0(find.package('gotchaslim'), '/gotcha_labeling.py'), quiet = T
  )

  cli::cli_alert_success("Environment setup looks alright, nice one!")
}
