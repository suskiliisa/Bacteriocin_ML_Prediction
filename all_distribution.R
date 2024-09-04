### FINAL OPTIMIZED VERSION #########

library(Biostrings)
library(seqinr)
library(protr)
library(stringr)
library(data.table)
library(parallel)

num_cores <- 16  

# Function to remove stop codons
remove_stop_codon <- function(seq) {
  return(gsub("\\*", "", seq))
}

# Function to extract distribution features
extractCTDD_revised <- function(x) {
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
  
  # Compute Distribution
  D <- vector("list", 8)
  for (i in 1:8) D[[i]] <- matrix(ncol = 5, nrow = 3)
  
  for (i in 1:8) {
    for (j in 1:3) {
      inds <- which(G[[i]] == paste0("G", j))
      quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
      
      quartiles[which(quartiles <= 0)] <- 1
      
      D[[i]][j, ] <- if (length(inds) > 0) {
        (inds[c(1, quartiles, length(inds))]) * 100 / n
      } else {
        0
      }
    }
  }
  
  D <- do.call(rbind, D)
  D <- as.vector(t(D))
  
  names(D) <- paste(
    rep(paste("prop", 1:8, sep = ""), each = 15),
    rep(rep(c(".G1", ".G2", ".G3"), each = 5), times = 8),
    rep(paste(".residue", c("0", "25", "50", "75", "100"), sep = ""), times = 24),
    sep = ""
  )
  
  flag <- matrix(0, 120, 1)
  lc <- c()
  
  for (i in seq(1, 120, 15)) {
    for (j in seq(i, i + 14)) {
      if (flag[j, 1] == 0) {
        z <- j
        count <- 1
        while (count <= 3) {
          lc <- c(lc, D[[z]])
          flag[z, 1] <- 1
          z <- z + 5
          count <- count + 1
        }
      }
    }
  }
  
  # Return the first 105 features
  return(lc[1:105])
}

# Function to process sequences in parallel
process_sequences <- function(index, ncrna2) {
  seq_id <- names(ncrna2)[index]
  x <- ncrna2[[index]]
  x <- remove_stop_codon(x)  # Remove stop codon
  d <- extractCTDD_revised(x)
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
  col_names <- c("ID", paste0("dist_", 1:105))  # 105 features
  setnames(results_dt, col_names)
  
  # Write to files
  output_prefix <- basename(dir)
  fwrite(results_dt, file.path("/nobackup/kbxz52/R_outputs/feat_extraction/DISTRIBUTION", paste0(output_prefix, "_distribution.txt")), sep = " ", col.names = TRUE)
  fwrite(results_dt, file.path("/nobackup/kbxz52/R_outputs/feat_extraction/DISTRIBUTION", paste0(output_prefix, "_distribution.csv")))
} else {
  warning("File not found: ", file_path)
}

