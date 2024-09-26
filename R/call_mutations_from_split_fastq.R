call_mutations_from_split_fastq <- function(
    filtered_files_in = 'path/to/filtered/fastq/',
    single_cell_csv = 'path_to/singlecell.csv',
    wt_max_mismatch = 0,
    mut_max_mismatch = 0,
    which_read = "R1",
    primer_sequence,
    primer_max_mismatch = 3,
    wt_sequence,
    mut_sequence,
    mutation_start,
    mutation_end,
    batch_size = 20,                     # Number of chunks per batch
    output_dir = 'path/to/output_batches/', # Directory to save batch outputs
    checkpoint = FALSE                    # Whether to use checkpoints
){
  
  cli::cli_h1("Starting mutation calling from split FASTQ files in batches")
  
  # Check if barcodes file exists
  if (!file.exists(single_cell_csv)) {
    cli::cli_abort("ATAC barcodes are not accessible under {single_cell_csv}")
  }
  
  if (batch_size %% 3 != 0) {
    cli::cli_alert_warning("Batch size is not a multiple of 3. Adjusting to nearest multiple.")
    batch_size <- ceiling(batch_size / 3) * 3  # Adjust to the nearest multiple of 3
  }
  
  # List all FASTQ files with full paths
  fastq_files <- list.files(path = filtered_files_in, pattern = ".gz", full.names = TRUE)
  
  # Check if any FASTQ files are found
  if (length(fastq_files) == 0) {
    cli::cli_abort("No .fastq.gz files found under {filtered_files_in}")
  }
  
  cli::cli_alert_info("Found {length(fastq_files)} .fastq.gz files in {filtered_files_in}")
  
  # Extract chunk number (e.g., "chunk01") and R1/R2/R3 label
  chunk_info <- gsub(".*(chunk[0-9]+).*", "\\1", basename(fastq_files))
  read_info <- gsub(".*(R[1-3]).*", "\\1", basename(fastq_files))
  
  # Combine into a dataframe for easier manipulation
  fastq_df <- data.frame(
    file = fastq_files,
    chunk = chunk_info,
    read = read_info,
    stringsAsFactors = FALSE
  )
  
  # Ensure that for each chunk, we have R1, R2, and R3
  valid_chunks <- fastq_df %>%
    dplyr::group_by(chunk) %>%
    dplyr::filter(sum(read == "R1") > 0 & sum(read == "R2") > 0 & sum(read == "R3") > 0) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(chunk, read) 
  
  # Now split these valid chunks into batches of the desired size
  batches <- split(valid_chunks$file, ceiling(seq_along(valid_chunks$file) / batch_size))
  
  total_batches <- length(batches)
  
  # Handle checkpoint logic
  if (checkpoint) {
    # List already existing CSV files in the output directory
    batch_output_files <- list.files(path = output_dir, pattern = "\\.csv$", full.names = TRUE)
    
    if (length(batch_output_files) > 0) {
      cli::cli_alert_info("Checkpoint enabled: Found existing batch output files.")
    } else {
      cli::cli_alert_info("No existing batch files found, proceeding with new processing.")
    }
  } else {
    # Remove existing output directory and create a new one
    if (dir.exists(output_dir)) {
      unlink(output_dir, recursive = TRUE) # Remove existing files_out directory
    }
    dir.create(output_dir, recursive = TRUE) # Create a new output directory
    cli::cli_alert_info("Created fresh output directory at {output_dir}")
  }
  
  cli::cli_alert_info("Determining cell barcodes to genotype")
  
  barcodes <- read.csv(single_cell_csv, header = TRUE, stringsAsFactors = FALSE)
  barcodes <- barcodes %>% dplyr::filter(barcode != 'NO_BARCODE') %>% dplyr::filter(excluded_reason == 0)
  barcodes$barcode <- substr(barcodes$barcode, 1, 16)
  barcodes <- barcodes[-1, 1]
  
  rev_barcodes <- reverse_complement(barcodes)
  
  batch_output_files <- vector("list", total_batches)
  
  # Iterate over each batch
  for(batch_idx in seq_along(batches)){
    # Generate output file path for the batch
    batch_output_file <- file.path(output_dir, sprintf("merged_output_batch_%03d.csv", batch_idx))
    
    if (checkpoint && file.exists(batch_output_file)) {
      cli::cli_alert_info("Skipping batch {batch_idx}/{total_batches}, checkpoint file exists.")
      batch_output_files[[batch_idx]] <- batch_output_file
      next
    }
    
    cli::cli_alert_info("Processing batch {batch_idx}/{total_batches}")
    current_batch <- batches[[batch_idx]]
    
    # Read FASTQ files
    progressr::with_progress({
      p <- progressr::progressor(along = current_batch)
      fastq_files_sequences <- future_lapply(current_batch, function(file) {
        p()
        temp <- ShortRead::readFastq(file)
        data.frame(
          Identifier = sub(" .*", "", as.character(ShortRead::id(temp))),
          Sequence = as.character(ShortRead::sread(temp)),
          stringsAsFactors = FALSE
        )
      })
    })
    
    names(fastq_files_sequences) <- basename(current_batch)
    names(fastq_files_sequences) <- sub("(_filtered\\.fastq\\.gz|\\.fastq\\.gz)$", "", names(fastq_files_sequences))
    
    filter_by_chunk <- function(chunk_index, barcodes) {
      R1_chunk <- as.data.table(fastq_files_sequences[grepl("R1", names(fastq_files_sequences))][[chunk_index]])
      R2_chunk <- as.data.table(fastq_files_sequences[grepl("R2", names(fastq_files_sequences))][[chunk_index]])
      R3_chunk <- as.data.table(fastq_files_sequences[grepl("R3", names(fastq_files_sequences))][[chunk_index]])
      
      setkey(R1_chunk, Identifier)
      setkey(R2_chunk, Identifier)
      setkey(R3_chunk, Identifier)
      
      R1_filtered <- R1_chunk[R3_chunk, nomatch = 0]  # Join R2 with R1
      R2_filtered <- R2_chunk[R3_chunk, nomatch = 0]  # Join R3 with R1
      
      combined_df <- data.table(
        R1 = R1_filtered$Sequence,
        R2 = R2_filtered$Sequence,
        R3 = R3_chunk$Sequence
      )
      
      combined_df_filtered <- combined_df[R2 %in% rev_barcodes]
      
      return(combined_df_filtered)
    }
    
    r1_files <- names(fastq_files_sequences[grepl("R1", names(fastq_files_sequences))])
    
    length <- length(r1_files)    
    
    progressr::with_progress({
      p <- progressr::progressor(along = seq_along(r1_files))
      filtered_reads <- future_lapply(
        seq_along(r1_files), 
        function(i) {
          p()
          filter_by_chunk(i)
        }
      )
    })
    
    names(filtered_reads) <- sprintf("chunk%02d", seq_along(filtered_reads))
    
    rm(fastq_files_sequences)
    
    gc()
    
    # Filter reads that match the primer sequence and apply reverse complement
    progressr::with_progress({
      p <- progressr::progressor(along = filtered_reads)
      filtered_reads <- future_lapply(seq_along(filtered_reads), function(i) {
        p()
        x <- filtered_reads[[i]]
        
        # Filter reads based on primer matching
        matches <- vcountPattern(primer_sequence, x[[which_read]], 
                                 max.mismatch = primer_max_mismatch, 
                                 with.indels = FALSE) >= 1
        
        # Apply filtering and reverse complement to the sequences
        filtered_reads <- x[matches, ]

        return(filtered_reads)
      })
    })
    
    cli::cli_alert_info("Counting wildtype and mutant reads in batch {batch_idx}")
    
    # Count wildtype and mutant reads
    progressr::with_progress({
      p <- progressr::progressor(along = filtered_reads)
      genotypes <- future_lapply(filtered_reads, function(x) {
        # Extract the relevant portion of the read for comparison
        read_subseq <- substr(as.character(x[[which_read]]), start = mutation_start, stop = mutation_end)
        x$wildtype <- vcountPattern(wt_sequence, read_subseq, max.mismatch = wt_max_mismatch, with.indels = FALSE) >=1
        x$mutant <- vcountPattern(mut_sequence, read_subseq, max.mismatch = mut_max_mismatch, with.indels = FALSE) >=1
        
        # Initialize genotype column
        x$ReadGenotype <- NA
        
        # Assign genotypes based on patterns
        x$ReadGenotype[x$wildtype & x$mutant] <- "ambiguous"
        x$ReadGenotype[!x$wildtype & !x$mutant] <- "ambiguous"
        x$ReadGenotype[x$wildtype & !x$mutant] <- "wildtype"
        x$ReadGenotype[!x$wildtype & x$mutant] <- "mutant"
        
        return(x)
      }, future.globals = list(
        wt_sequence = wt_sequence, 
        mut_sequence = mut_sequence, 
        which_read = which_read, 
        mutation_start = mutation_start, 
        mutation_end = mutation_end, 
        wt_max_mismatch = wt_max_mismatch, 
        mut_max_mismatch = mut_max_mismatch
      ), 
      future.packages = c("Biostrings"))
    })

    names(genotypes) <-  sprintf("chunk%02d", seq_along(genotypes))
    
    progressr::with_progress({
      p <- progressr::progressor(along = genotypes)
      counts_per_chunk <- future_lapply(genotypes, function(x) {
        p()
        x <- x[, c("R2", "wildtype", "mutant", "ReadGenotype")]
        
        # Summarise counts
        x_summary <- x %>%
          group_by(R2) %>%
          summarise(
            wildtype_count = sum(ReadGenotype == "wildtype"),
            mutant_count = sum(ReadGenotype == "mutant"),
            ambiguous_count = sum(ReadGenotype == "ambiguous"),
            .groups = 'drop'  # Avoid warnings about grouping
          )
        
        return(x_summary)  # Return as a data frame
      })
    })
    
    counts <- dplyr::bind_rows(counts_per_chunk)  # Combine data frames
    
    # Aggregate counts
    counts <- counts %>%
      dplyr::group_by(R2) %>%
      dplyr::summarise(
        wildtype_count = sum(wildtype_count, na.rm = TRUE),
        mutant_count = sum(mutant_count, na.rm = TRUE),
        ambiguous_count = sum(ambiguous_count, na.rm = TRUE),
        .groups = 'drop'  # Avoid warnings about grouping
      ) %>%
      tidyr::replace_na(list(
        wildtype_count = 0,
        mutant_count = 0,
        ambiguous_count = 0
      )) 
    
    colnames(counts) <- c('barcode', 'wildtype_count', 'mutant_count', 'ambiguous_count')
    
    # Save the merged output for the current batch
    write.csv(counts, batch_output_file, row.names = FALSE)
    batch_output_files[[batch_idx]] <- batch_output_file
    
    # Clean up memory
    gc()
    
    cli::cli_alert_success("Completed processing batch {batch_idx}")
  }
  
  cli::cli_alert_success("Finished genotyping per batch!")
}
