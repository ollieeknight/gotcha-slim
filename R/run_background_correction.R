run_background_correction <- function(input_file, output_dir, saturation = FALSE) {
  # Source the Python script
  py_run_file("~/work/scripts/gotcha-slim/python/gotcha_labeling.py")
  
  # Call the GotchaLabeling function from Python
  result <- py$GotchaLabeling(merged_genotype_file = input_file,
                              saturation = saturation)
  
  # Save the result to the specified output directory
  output_file <- file.path(output_dir, "genotype_labels.csv")
  result$to_csv(output_file, index = FALSE)
  
  return(result)
}
