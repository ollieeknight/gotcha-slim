reverse_complement <- function(sequences) {
  dna_strings <- Biostrings::DNAStringSet(sequences)
  rev_complements <- reverseComplement(dna_strings)
  return(as.character(rev_complements))
}