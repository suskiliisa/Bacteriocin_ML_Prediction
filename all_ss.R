### FINAL OPTIMIZED VERSION #########

library(seqinr)
library(Biostrings)
library(DECIPHER)
library(stringr)
library(data.table)
library(parallel)

# Set the number of cores to use
num_cores <- 16  # Adjust this number based on your supercomputer's available resources

# Function to remove stop codons
remove_stop_codon <- function(seq) {
  return(gsub("\\*", "", seq))
}

# Process sequences in parallel and extract secondary structure features
process_sequences <- function(index, ncrna2, hec) {
  seq_id <- names(ncrna2)[index]
  seq <- ncrna2[[index]]
  seq <- remove_stop_codon(seq)
  
  t <- hec[index]
  data_psi <- strsplit(t, "")[[1]]
  len <- length(data_psi)
  add_ch <- ""
  add_ch_exclude <- ""
  flag <- ""
  freq <- 0
  SH <- 0
  SE <- 0
  SC <- 0
  
  for (v in 1:len) {
    if (toString(data_psi[v]) == "H") {
      SH <- SH + v
    }
    if (toString(data_psi[v]) == "E") {
      SE <- SE + v
    }
    if (toString(data_psi[v]) == "C") {
      SC <- SC + v
    }
    add_ch <- paste0(add_ch, toString(data_psi[v]))
  }
  
  rr <- rle(strsplit(add_ch, "")[[1]])
  MH <- max(rr$lengths[which(rr$values == "H")], na.rm = TRUE)
  if (is.infinite(MH)) MH <- 0
  CMVH <- SH / (len * (len - 1))
  NMH <- MH / len
  
  ME <- max(rr$lengths[which(rr$values == "E")], na.rm = TRUE)
  if (is.infinite(ME)) ME <- 0
  CMVE <- SE / (len * (len - 1))
  NME <- ME / len
  
  MC <- max(rr$lengths[which(rr$values == "C")], na.rm = TRUE)
  CMVC <- SC / (len * (len - 1))
  
  rr_len <- length(rr$values)
  for (v in 1:rr_len) {
    if (rr$values[v] != "C") {
      if (rr$values[v] != flag) {
        add_ch_exclude <- paste0(add_ch_exclude, toString(rr$values[v]))
        flag <- toString(rr$values[v])
        freq <- freq + 1
      }
    }
  }
  
  count_EHE <- gregexpr("(?=EHE)", add_ch_exclude, perl = TRUE)
  count_EHE1 <- count_EHE[[1]]
  if (count_EHE1[1] < 0) {
    count_EHE_len <- 0
  } else {
    count_EHE_len <- length(count_EHE1)
  }
  
  if (freq > 2) {
    f_EHE <- count_EHE_len / (freq - 2)
  } else {
    f_EHE <- 0
  }
  
  temp <- c(seq_id, CMVH, CMVE, CMVC, NMH, NME, f_EHE)
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
  
  # Load the secondary structure prediction
  aa <- readAAStringSet(file_path)
  hec <- PredictHEC(aa)
  
  # Process sequences in parallel
  results <- mclapply(1:length(ncrna2), process_sequences, ncrna2 = ncrna2, hec = hec, mc.cores = num_cores)
  
  # Remove NULL entries from results
  results <- do.call(rbind, results[!sapply(results, is.null)])
  
  # Convert to data.table for efficient writing
  results_dt <- as.data.table(results)
  
  # Constructing column names
  col_n <- "ss_1"
  for (k in 2:6) {
    col_n <- paste0(col_n, " ", "ss_", k)
  }
  col_names <- c("ID", unlist(strsplit(col_n, " ")))
  
  setnames(results_dt, old = names(results_dt), new = col_names)
  
  # Define output directories and ensure they exist
  output_dir <- "/nobackup/kbxz52/R_outputs/feat_extraction/SS"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write to files
  output_prefix <- basename(dir)
  fwrite(results_dt, file.path(output_dir, paste0(output_prefix, "_ss_features.txt")), sep = " ", col.names = TRUE)
  fwrite(results_dt, file.path(output_dir, paste0(output_prefix, "_ss_features.csv")))
} else {
  warning("File not found: ", file_path)
}

