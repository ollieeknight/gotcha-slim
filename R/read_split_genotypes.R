read_split_genotypes <- function(files_in, existing_merged_file = TRUE) {
  
  merged_file <- file.path(files_in, "merged_genotypes.csv")
  
  if (existing_merged_file && file.exists(merged_file)) {
    message("Reading from existing file")
    
    called_genotypes <- read.csv(merged_file, stringsAsFactors = FALSE)
    
  } else {
    message("Generating a new merged file")
    
    files <- list.files(files_in, pattern = "batch_.*\\.csv$", full.names = TRUE)
    
    read_genotype_file <- function(file) {
      df <- read.csv(file, stringsAsFactors = FALSE)

            return(df)
    }
    
    all_genotypes <- do.call(rbind, lapply(files, read_genotype_file))
    
    called_genotypes <- all_genotypes %>%
      dplyr::group_by(barcode) %>%
      dplyr::summarise(
        wildtype_count = sum(wildtype_count, na.rm = TRUE),
        mutant_count = sum(mutant_count, na.rm = TRUE),
        ambiguous_count = sum(ambiguous_count, na.rm = TRUE)
        
      ) %>%
      dplyr::ungroup()
  
    final_genotypes <- called_genotypes %>%
      dplyr::select(barcode, wildtype_count, mutant_count, ambiguous_count)
    
    # Write the merged genotypes to a CSV file
    write.csv(final_genotypes, file = merged_file, row.names = FALSE)
    message(paste("Merged genotypes saved to", merged_file))
    
    # Return the merged dataframe
    return(final_genotypes)
  }
}
