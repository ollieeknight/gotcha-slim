Sys.setenv(PATH = paste('~/work/bin/miniconda3/envs/r-reticulate/lib/python3.10/site-packages/', Sys.getenv()['PATH'], sep = ':'))
library(reticulate)
assignInNamespace('is_conda_python', function(x){ return(FALSE) }, ns = 'reticulate')
use_miniconda('~/work/bin/miniconda3/envs/r-reticulate/', required = T)
set.seed(123)

library(future)
library(future.apply)
library(ShortRead)
library(fs)
library(cli)
library(future)
library(stringr)
library(progressr)
library(data.table)
library(Biostrings)
library(dplyr)
library(ggplot2)

options(future.globals.maxSize = Inf)  # Adjust size limit

plan(multisession, workers = 32)

scripts <- list.files('~/work/scripts/gotcha-slim/R/', pattern = "\\.R$", full.names = TRUE)

lapply(scripts, source)

check_gotcha_slim_environment()

merge_lane_fastq_files(files_in = '~/share/ngs/raw/ASAP_SI_GENO/outs/fastq_path/',
                        files_out = '~/share/ngs/raw/ASAP_SI_GENO/01merged/')

split_merged_fastq_files(files_in = '~/share/ngs/raw/ASAP_SI_GENO/01merged/',
                  files_out = '~/share/ngs/raw/ASAP_SI_GENO/02split/', 
                  reads = 1000000)

filter_split_fastq_files(files_in = '~/share/ngs/raw/ASAP_SI_GENO/02split/',
                         files_out = '~/share/ngs/raw/ASAP_SI_GENO/03filtered/', 
                         min_quality = 15, 
                         min_bases = 1, 
                         which_read = 'R3',
                         read_region = c(69:71)
                   )

call_mutations_from_split_fastq(
  filtered_files_in = '~/share/ngs/raw/ASAP_SI_GENO/03filtered/',
  single_cell_csv = '~/work/data/adaptive/data/ASAP_GoTChA_HC13_exp1_libA/outs/singlecell.csv',
  output_dir = '~/share/ngs/raw/ASAP_SI_GENO/04genotype/',
  checkpoint = F,
  batch_size = 12,
  wt_max_mismatch = 0,
  mut_max_mismatch = 0,
  which_read = "R3",
  primer_sequence = 'GGCAGCCCTATGATTGGAGTAG',
  primer_max_mismatch = 3,
  wt_sequence =  "ACG",
  mut_sequence = "ATG",
  mutation_start = 69, 
  mutation_end = 71
)

called_genotypes <- read_split_genotypes(
  files_in = '~/share/ngs/raw/ASAP_SI_GENO/04genotype/', 
  existing_merged_file = F
)

df <- called_genotypes
df$barcode <- reverse_complement(df$barcode)

real <- readLines('adaptive/data/ASAP_GoTChA_HC13_exp1_libA/outs/filtered_peak_bc_matrix/barcodes.tsv')
real <- substr(real, 1, 16)

df <- df %>%
  mutate(real = barcode %in% real)

clr_function <- function(x) {
  log1p(x / (exp(sum(log1p(x[x > 0]), na.rm = TRUE) / length(x))))
}

# Apply the CLR-like normalisation to the 'values' column and add it as a new column
called_genotypes$wildtype_clr <- clr_function(called_genotypes$wildtype_count)
called_genotypes$mutant_clr <- clr_function(called_genotypes$mutant_count)

barcodes <- readLines('~/work/data/adaptive/data/ASAP_GoTChA_HC13_exp1_libA/outs/filtered_peak_bc_matrix/barcodes.tsv')
barcodes <- gsub('-1', '', barcodes)

called_genotypes$barcode <- reverse_complement(called_genotypes$barcode)

called_genotypes <- called_genotypes[called_genotypes$barcode %in% barcodes,]

ggplot(called_genotypes, aes(x = wildtype_clr, y = mutant_clr)) +
  geom_point() +
  labs(
    title = NULL,
    x = "Wildtype CLR",
    y = "Mutant CLR"
  ) +
  theme_classic() + coord_fixed() +
  geom_hline(yintercept = 2.5) +
  geom_vline(xintercept = 2.5)
