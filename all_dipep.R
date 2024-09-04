### FINAL OPTIMIZED VERSION #########

library(Biostrings)
library(seqinr)
library(protr)
library(data.table)
library(parallel)

# Set the number of cores to use
num_cores <- 16  # Adjust this number based on your supercomputer's available resources

# Function to remove stop codons
remove_stop_codon <- function(seq) {
  return(gsub("\\*", "", seq))
}

# Dipeptide Composition features extraction function
extractDC <- function(x) {
  x <- toupper(x)
  aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  dipeptides <- outer(aa, aa, paste0)
  dipep_dict <- setNames(rep(0, length(dipeptides)), dipeptides)
  
  for (i in 1:(nchar(x) - 1)) {
    dipep <- substr(x, i, i + 1)
    if (dipep %in% names(dipep_dict)) {
      dipep_dict[dipep] <- dipep_dict[dipep] + 1
    }
  }
  
  dipep_freq <- dipep_dict / (nchar(x) - 1)
  return(dipep_freq)
}

# Function to process sequences in parallel
process_sequences <- function(index, ncrna2) {
  seq_id <- names(ncrna2)[index]
  x <- ncrna2[[index]]
  x <- remove_stop_codon(x)
  d <- extractDC(x)
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

# Importing reads
file_path <- file.path(dir, "predicted_proteins.faa")
if (file.exists(file_path)) {
  ncrna2 <- read.fasta(file = file_path, as.string = TRUE, seqtype = "AA")
  
  # Process sequences in parallel
  results <- mclapply(1:length(ncrna2), process_sequences, ncrna2 = ncrna2, mc.cores = num_cores)
  
  # Remove NULL entries from results
  results <- do.call(rbind, results[!sapply(results, is.null)])
  
  # Convert to data.table for efficient writing
  results_dt <- as.data.table(results)
  
  # Constructing column names
  col_n <- "dipep_1"
  for (k in 2:400) {
    col_n <- paste0(col_n, " ", "dipep_", k)
  }
  col_names <- c("ID", unlist(strsplit(col_n, " ")))
  
  setnames(results_dt, old = names(results_dt), new = col_names)
  
  # Define output directories and ensure they exist
  output_dir <- "/nobackup/kbxz52/R_outputs/feat_extraction/DIPEP"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write to files
  output_prefix <- basename(dir)
  fwrite(results_dt, file.path(output_dir, paste0(output_prefix, "_dipep_features.txt")), sep = " ", col.names = TRUE)
  fwrite(results_dt, file.path(output_dir, paste0(output_prefix, "_dipep_features.csv")))
} else {
  warning("File not found: ", file_path)
}

