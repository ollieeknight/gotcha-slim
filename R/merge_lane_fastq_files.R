merge_lane_fastq_files <- function(files_in, files_out) {
  if (!dir_exists(files_in)) stop("The specified input path does not exist.")
  if (!dir_exists(files_out)) dir_create(files_out)
  
  # List all fastq.gz files in the input folder
  file.index <- dir_ls(files_in, regexp = "\\.fastq\\.gz$")
  
  if (length(file.index) == 0) {
    cli::cli_abort("No .fastq.gz files found under {.path {files_in}}.")
  }
  
  # Extract prefixes before lane and read type
  prefixes <- str_extract(path_file(file.index), "^.*(?=_S[0-9]{1,2}_L[0-9]{3}_R[1-3])")
  
  grouped_files <- split(file.index, list(prefixes, str_extract(path_file(file.index), "_R[1-3]")))
  
  # Clean the names by removing the dot
  names(grouped_files) <- gsub("\\.", "", names(grouped_files))
  
  # Function to merge files with the same prefix and read type
  merge_fastq <- function(files, prefix) {
    clean_prefix <- str_remove(prefix, "\\.$")  # Remove any trailing dot
    merged_file <- file.path(files_out, paste0(clean_prefix, "_merged.fastq.gz"))
    
    system(paste("cat", paste(files, collapse = " "), ">", merged_file))
    cli::cli_alert_success("{merged_file} created successfully.")
  }
  
  # Parallelised merging with future_lapply
  invisible(future_lapply(names(grouped_files), function(name) {
    merge_fastq(grouped_files[[name]], name)
  }))
  
  cli::cli_alert_success("All fastq.gz files merged")
}
