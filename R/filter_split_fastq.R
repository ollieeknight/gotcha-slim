filter_split_fastq_files = function(files_in,
                                    files_out,
                                    min_quality = 15,
                                    min_bases = 1,
                                    which_read = "R1",
                                    read_region = NULL) {
  
  # List all FASTQ files
  fastq_files = dir(path = files_in, pattern = ".gz", full.names = TRUE)
  
  # Check and inform based on the number of files found
  if (length(fastq_files) == 0) {
    cli::cli_abort("No .fastq.gz files found under {.path {files_in}}.")
  } else {  
    cli::cli_h1("Starting to filter split .fastq.gz files")
    cli::cli_alert_info("Found {length(fastq_files)} .fastq.gz file(s) under {.path {files_in}}.")
  }
  
  # Ensure output directory exists
  fs::dir_create(files_out)
  
  processed_files <- character(0)
  
  # Read FASTQ files in the current batch
  progressr::with_progress({
    p <- progressr::progressor(along = fastq_files)
    results <- future_lapply(fastq_files, function(file) {
      # Filter based on which read
      if (!grepl(which_read, basename(file))) {
        p()
        return(NULL)
      }
      
      # Read and filter quality
      temp = ShortRead::readFastq(file)
      q = as(quality(temp), "matrix")
      if (!is.null(read_region) && min_bases <= length(read_region)) {
        q = q[, read_region]
      }
      
      index = rowSums(q <= min_quality) < min_bases
      filtered_data = temp[index]
      
      output_file = file.path(files_out, 
                              sub("(.+)(_chunk\\d+)(\\.fastq\\.gz)$", "\\1\\2_filtered\\3", basename(file)))
      
      ShortRead::writeFastq(object = filtered_data, file = output_file, mode = "w")
      
      processed_files <<- c(processed_files, file)
      p()
    })
  })
  
  unprocessed_files <- setdiff(fastq_files, processed_files)
  number = length(unprocessed_files)
  total = length(fastq_files)
  
  if (number > 0) {
    cli::cli_alert_info("{number} files were not processed")
    cli::cli_alert_info("Copying unprocessed files to output directory...")
    
    progressr::with_progress({
      p <- progressr::progressor(along = unprocessed_files)
      future_lapply(unprocessed_files, function(file) {
        if (!grepl(which_read, basename(file))) {
          fs::file_copy(file, file.path(files_out, basename(file)), overwrite = FALSE)
        }
        p()
      })
    })
    
    cli::cli_alert_success("Processed {total - number} out of {total} files")
  }
}
