### FINAL OPTIMIZED VERSION #########

library(Biostrings)
library(seqinr)
library(protr)
library(data.table)
library(parallel)

# Set the number of cores to use
num_cores <- 10  

# Function to remove stop codons
remove_stop_codon <- function(seq) {
  return(gsub("\\*", "", seq))
}

# Modified function to store problematic sequences
extractAAC_BAC <- function(x) {
  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )
  
  if (any(!strsplit(x, split = "")[[1]] %in% AADict)) {
    warning("Sequence contains unrecognized amino acid type")
    return(NULL)
  }
  
  AAC <- summary(
    factor(strsplit(x, split = "")[[1]], levels = AADict),
    maxsum = 21
  ) / nchar(x)
  
  return(AAC)
}

# Function to process sequences in parallel
process_sequences <- function(index, ncrna2) {
  seq_id <- names(ncrna2)[index]
  x <- ncrna2[[index]]
  x <- remove_stop_codon(x)  # Remove stop codon
  d <- extractAAC_BAC(x)
  if (is.null(d)) {
    return(NULL)
  }
  temp <- c(seq_id, as.numeric(d))
  return(temp)
}

# Get the directory from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No directory provided")
}
dir <- args[1]

file_path <- file.path(dir, "predicted_proteins.faa")
if (file.exists(file_path)) {
  ncrna2 <- read.fasta(file = file_path, as.string = TRUE, seqtype = "AA")
  results <- mclapply(1:length(ncrna2), process_sequences, ncrna2, mc.cores = num_cores)
  results <- do.call(rbind, results[!sapply(results, is.null)])
  results_dt <- as.data.table(results)
  col_names <- c("ID", paste0("aac_", 1:20))
  setnames(results_dt, col_names)
  
  # Write to files
  output_prefix <- basename(dir)
  fwrite(results_dt, file.path("/nobackup/kbxz52/R_outputs/feat_extraction", paste0(output_prefix, "_acc.txt")), sep = " ", col.names = TRUE)
  fwrite(results_dt, file.path("/nobackup/kbxz52/R_outputs/feat_extraction", paste0(output_prefix, "_acc.csv")))
} else {
  warning("File not found: ", file_path)
}

