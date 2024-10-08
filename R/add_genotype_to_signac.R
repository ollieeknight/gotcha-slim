add_genotypes_to_signac_object <- function(alldata,
                                called_genotypes,
                                gene_id = '') {
  
  if (!inherits(alldata, "Seurat")) {
    stop("The specified alldata object is not a Seurat class")
  }
  if (length(samples) < length(gotchaList)) {
    stop("Provide names for all samples in the list of GoTChA object paths")
  }
  if (length(samples) > length(gotchaList)) {
    stop("Provide names only for those samples present in the GoTChA object paths, in the correct order")
  }
  
  message("------- LOADING GOTCHA OUTPUTS -------")
  geno <- lapply(gotchaList, function(x) {
    return(get(load(x)))
  })
  names(geno) <- samples
  
  message("------- ADDING SAMPLE AND -1 TO THE BARCODES TO MATCH SIGNAC FORMAT ------")
  geno <- lapply(samples, function(x) {
    geno[[x]]$WhiteListMatch <- paste0(x, "#", geno[[x]]$WhiteListMatch, "-1")
    rownames(geno[[x]]) <- geno[[x]]$WhiteListMatch
    return(geno[[x]])
  })
  names(geno) <- samples
  
  # Check how many barcodes in genotyping library match the ones in the Signac object
  detected <- sum(unlist(lapply(geno, function(x) x$WhiteListMatch %in% colnames(alldata))))
  if (sum(detected == 0)) {
    stop("No matching barcodes detected! Check for reverse complement in MutationCalling function")
  }
  message(paste0(detected, " total matching barcodes detected in Signac object (", 
                 round(detected / ncol(alldata) * 100, 2), "%)"))
  
  message("------- PERFORMING SAMPLE-AWARE NOISE CORRECTION BASED ON EMPTY DROPLETS -------")
  
  geno.corrected <- lapply(samples, function(i) {
    message(paste0("Processing sample: ", i))
    
    empty.droplets <- geno[[i]][!(geno[[i]]$WhiteListMatch %in% colnames(alldata)), ]
    message(paste0(nrow(empty.droplets), " empty droplets will be used for noise correction"))
    
    noise.wt <- quantile(empty.droplets$WTcount, probs = seq(0, 1, 0.01))[percentile.use + 1]
    noise.mut <- quantile(empty.droplets$MUTcount, probs = seq(0, 1, 0.01))[percentile.use + 1]
    
    if (plot) {
      hist(log10(empty.droplets$WTcount + 1), breaks = 100, 
           main = paste0(i, " - WT reads in empty droplets"), xlab = "Wildtype read count (log10)")
      abline(v = log10(noise.wt + 1), lty = 2)
      hist(log10(empty.droplets$MUTcount + 1), breaks = 100, 
           main = paste0(i, " - MUT reads in empty droplets"), xlab = "Mutant read count (log10)")
      abline(v = log10(noise.mut + 1), lty = 2)
    }
    
    message(paste0("Noise in WT reads: ", round(noise.wt, 4), 
                   " --- Noise in MUT reads: ", round(noise.mut, 4), 
                   " --- Quantile used: ", names(noise.wt)))
    
    geno[[i]]$WTcount_noise_corrected <- geno[[i]]$WTcount - noise.wt
    geno[[i]]$WTcount_noise_corrected[geno[[i]]$WTcount_noise_corrected < 0] <- 0
    geno[[i]]$MUTcount_noise_corrected <- geno[[i]]$MUTcount - noise.mut
    geno[[i]]$MUTcount_noise_corrected[geno[[i]]$MUTcount_noise_corrected < 0] <- 0
    
    return(geno[[i]])
  })
  
  geno <- as.data.frame(do.call(rbind, geno.corrected))
  rownames(geno) <- geno$WhiteListMatch
  
  message("------- ADDING GENOTYPING INFORMATION TO SIGNAC OBJECT METADATA ------")
  
  message("Subsetting genotyping data to cells present in the Signac object...")
  geno.sub <- geno[colnames(alldata), ]
  
  message("Adding raw GoTChA read counts...")
  alldata$meta.data[, paste0(target, "_WTcount")] <- geno.sub[colnames(alldata), ]$WTcount
  alldata$meta.data[, paste0(target, "_MUTcount")] <- geno.sub[colnames(alldata), ]$MUTcount
  
  alldata$meta.data[, paste0(target, "_MUTfraction")] <- geno.sub[colnames(alldata), ]$MUTcount /
    (geno.sub[colnames(alldata), ]$MUTcount + geno.sub[colnames(alldata), ]$WTcount)
  alldata$meta.data[, paste0(target, "_TotalCounts")] <- geno.sub[colnames(alldata), ]$WTcount +
    geno.sub[colnames(alldata), ]$MUTcount
  
  message("Adding raw genotype calls...")
  alldata$meta.data[, paste0(target, "_Genotype")] <- NA
  alldata$meta.data[, paste0(target, "_Genotype")][alldata$meta.data[, paste0(target, "_MUTfraction")] <= genotyping.thresholds[1]] <- "WT"
  alldata$meta.data[, paste0(target, "_Genotype")][alldata$meta.data[, paste0(target, "_MUTfraction")] >= genotyping.thresholds[2]] <- "MUT"
  alldata$meta.data[, paste0(target, "_Genotype")][alldata$meta.data[, paste0(target, "_MUTfraction")] < genotyping.thresholds[2] & 
                                                     alldata$meta.data[, paste0(target, "_MUTfraction")] > genotyping.thresholds[1]] <- "HET"
  
  message("Adding noise-corrected GoTChA read counts...")
  alldata$meta.data[, paste0(target, "_WTcount_noise_corrected")] <- geno.sub[colnames(alldata), ]$WTcount_noise_corrected
  alldata$meta.data[, paste0(target, "_MUTcount_noise_corrected")] <- geno.sub[colnames(alldata), ]$MUTcount_noise_corrected
  
  alldata$meta.data[, paste0(target, "_MUTfraction_noise_corrected")] <- alldata$meta.data[, paste0(target, "_MUTcount_noise_corrected")] /
    (alldata$meta.data[, paste0(target, "_MUTcount_noise_corrected")] + alldata$meta.data[, paste0(target, "_WTcount_noise_corrected")])
  alldata$meta.data[, paste0(target, "_TotalCounts_noise_corrected")] <- alldata$meta.data[, paste0(target, "_MUTcount_noise_corrected")] +
    alldata$meta.data[, paste0(target, "_WTcount_noise_corrected")]
  
  message("Adding noise-corrected genotype calls...")
  alldata$meta.data[, paste0(target, "_Genotype_noise_corrected")] <- NA
  alldata$meta.data[, paste0(target, "_Genotype_noise_corrected")][alldata$meta.data[, paste0(target, "_MUTfraction_noise_corrected")] < genotyping.thresholds[1]] <- "WT"
  alldata$meta.data[, paste0(target, "_Genotype_noise_corrected")][alldata$meta.data[, paste0(target, "_MUTfraction_noise_corrected")] > genotyping.thresholds[2]] <- "MUT"
  alldata$meta.data[, paste0(target, "_Genotype_noise_corrected")][alldata$meta.data[, paste0(target, "_MUTfraction_noise_corrected")] < genotyping.thresholds[2] & 
                                                                     alldata$meta.data[, paste0(target, "_MUTfraction_noise_corrected")] > genotyping.thresholds[1]] <- "HET"
  
  message("------- DONE! ------")
  
  return(alldata)
}
