split_merged_fastq_files <- function(files_in, files_out, reads = 2000000) {
  
  # Get list of input FASTQ files
  fastq_files <- dir(files_in, pattern = ".gz", full.names = TRUE)
  
  if (length(fastq_files) == 0) {
    cli::cli_abort("No .fastq.gz files found under {.path {files_in}}.")
  } else {
    cli::cli_h1("Starting to split merged .fastq.gz files")
    cli::cli_inform("Found {length(fastq_files)} .fastq.gz file(s) under {.path {files_in}}.")    
  }
  
  # Create output folder if it doesn't exist
  if (!dir.exists(files_out)) {
    dir.create(files_out, recursive = TRUE)
  }
  # Read FASTQ files in the current batch
  progressr::with_progress({
    p <- progressr::progressor(along = fastq_files)
    future_lapply(fastq_files, function(x) {
      system(paste0(
        "gunzip -c ", x, 
        " | split - -l ", format(reads * 4, scientific = FALSE), 
        " --numeric-suffixes=1 --suffix-length=4 --filter='gzip -f - > $FILE.fastq.gz' ",
        files_out, gsub(".fastq.gz", "", basename(x)), "_chunk"
      ))
    })
  })
  
  cli::cli_alert_success("All .fastq.gz files successfully split")
}
