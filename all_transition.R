### FINAL OPTIMIZED VERSION #########

library(Biostrings)
library(seqinr)
library(protr)
library(data.table)
library(parallel)


num_cores <- 16 

# Function to remove stop codons
remove_stop_codon <- function(seq) {
  return(gsub("\\*", "", seq))
}

# Transition features extraction function
extractCTDT_revised <- function(x) {
  
  AADict <- c(
    "A", "C", "D", "E", "F", "G", "H", "I",
    "L", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" 
  )
  
  if (all(strsplit(x, split = "")[[1]] %in% AADict) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  group1 <- list(
    "hydrophobicity" = c("R", "K", "E", "D", "Q", "N"),
    "normwaalsvolume" = c("G", "A", "S", "T", "P", "D", "C"),
    "polarity" = c("L", "I", "F", "W", "C", "M", "V", "Y"),
    "polarizability" = c("G", "A", "S", "D", "T"),
    "charge" = c("K", "R"),
    "secondarystruct" = c("E", "A", "L", "M", "Q", "K", "R", "H"),
    "solventaccess" = c("A", "L", "F", "C", "G", "I", "V", "W"),
    "surfacetension" = c("G", "Q", "D", "N", "A", "H", "R")
  )
  
  group2 <- list(
    "hydrophobicity" = c("G", "A", "S", "T", "P", "H", "Y"),
    "normwaalsvolume" = c("N", "V", "E", "Q", "I", "L"),
    "polarity" = c("P", "A", "T", "G", "S"),
    "polarizability" = c("C", "P", "N", "V", "E", "Q", "I", "L"),
    "charge" = c(
      "A", "N", "C", "Q", "G", "H", "I", "L",
      "M", "F", "P", "S", "T", "W", "Y", "V"
    ),
    "secondarystruct" = c("V", "I", "Y", "C", "W", "F", "T"),
    "solventaccess" = c("R", "K", "Q", "E", "N", "D"),
    "surfacetension" = c("K", "T", "S", "E", "C")
  )
  
  group3 <- list(
    "hydrophobicity" = c("C", "L", "V", "I", "M", "F", "W"),
    "normwaalsvolume" = c("M", "H", "K", "F", "R", "Y", "W"),
    "polarity" = c("H", "Q", "R", "K", "N", "E", "D"),
    "polarizability" = c("K", "M", "H", "F", "R", "Y", "W"),
    "charge" = c("D", "E"),
    "secondarystruct" = c("G", "N", "P", "S", "D"),
    "solventaccess" = c("M", "S", "P", "T", "H", "Y"),
    "surfacetension" = c("I", "L", "M", "F", "P", "W", "Y", "V")
  )
  
  xSplitted <- strsplit(x, split = "")[[1]]
  n <- nchar(x)
  
  G <- vector("list", 8)
  for (i in 1:8) G[[i]] <- rep(NA, n)
  
  # Get groups for each property & each amino acid
  for (i in 1:8) {
    try(G[[i]][which(xSplitted %in% group1[[i]])] <- "G1")
    try(G[[i]][which(xSplitted %in% group2[[i]])] <- "G2")
    try(G[[i]][which(xSplitted %in% group3[[i]])] <- "G3")
  }
  
  # Combine single amino acids by a 2-length step
  for (i in 1:8) G[[i]] <- paste(G[[i]][-n], G[[i]][-1], sep = "")
  G <- lapply(G, function(x) factor(x, levels = c(
    "G1G2", "G2G1", "G1G3", "G3G1", "G2G3", "G3G2",
    "G1G1", "G2G2", "G3G3"
  )))
  
  GSummary <- lapply(G, summary)
  
  # Compute (n_rs + n_sr) / (N - 1)
  CTDT <- vector("list", 8)
  for (i in 1:8) {
    CTDT[[i]][1] <- sum(GSummary[[i]][c("G1G2", "G2G1")]) / (n - 1)
    CTDT[[i]][2] <- sum(GSummary[[i]][c("G1G3", "G3G1")]) / (n - 1)
    CTDT[[i]][3] <- sum(GSummary[[i]][c("G2G3", "G3G2")]) / (n - 1)
  }
  
  CTDT <- unlist(CTDT)
  
  return(CTDT[1:21])
}

# Function to process sequences in parallel
process_sequences <- function(index, ncrna2) {
  seq_id <- names(ncrna2)[index]
  x <- ncrna2[[index]]
  x <- remove_stop_codon(x)  # Remove stop codon
  d <- extractCTDT_revised(x)
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
  col_names <- c("ID", paste0("trans_", 1:21))  # 21 features
  setnames(results_dt, col_names)
  
  # Write to files
  output_prefix <- basename(dir)
  fwrite(results_dt, file.path("/nobackup/kbxz52/R_outputs/feat_extraction/TRANSITION", paste0(output_prefix, "_transition.txt")), sep = " ", col.names = TRUE)
  fwrite(results_dt, file.path("/nobackup/kbxz52/R_outputs/feat_extraction/TRANSITION", paste0(output_prefix, "_transition.csv")))
} else {
  warning("File not found: ", file_path)
}

