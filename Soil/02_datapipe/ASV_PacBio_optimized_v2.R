#!/usr/bin/env Rscript
# BiocManager::install(version = '3.20')
# Below are changes made to the original:

# 1. Removed primer removal step (checked and none of the fastq files contain the primers in the list)

# 2. Updated paths ..

# 3. For the error model, limit the data it processes
# err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=16, 
#                   multithread=TRUE, nbases=1e8)

# 4. Reduce BAND_SIZE from 32 to 16
# dd <- dada(drp, err=err, BAND_SIZE=16, multithread=TRUE)

# 5. Limit selfConsist iterations
# dada2::setDadaOpt(MAX_CONSIST=10)

# 6. Added timers for tracking progress throughout
# start_time <- Sys.time()
# err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=16, multithread=TRUE, nbases=1e8)
# end_time <- Sys.time()
# print(paste("Error learning completed in:", difftime(end_time, start_time, units="hours"), "hours"))


# ------------------------------------------------------------------------
# Preparation
# ------------------------------------------------------------------------
# Clear the object from memory
rm(list=ls())

# Start overall timer
script_start_time <- Sys.time()
print(paste("Script started at:", script_start_time))

#### CHANGE THIS PATHS ####
source <- "/lab/biohpc/students/sadiya/r_analysis_proj_2/project"
path.in <- "/lab/biohpc/students/sadiya/r_analysis_proj_2/project/01_demultiplex/out/" #needs a "/" at the end
path.primer <- "/lab/biohpc/students/sadiya/r_analysis_proj_2/project/01_demultiplex/" #needs a "/" at the end
path.rds <- "RDS/"
path.out <- "output/"
taxa_database <- "/lab/biohpc/students/sadiya/r_analysis_proj_2/project/02_datapipe/utax_reference_dataset_10.05.2021.fasta"

##########################

setwd(source)

# Libraries
print("Loading libraries...")
library_start_time <- Sys.time()
library(dada2)
library(Biostrings)
library(ShortRead)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(phyloseq)
library(xlsx)
library(seqinr)
library(vegan)
library(tibble)
library(dplyr)
library(DECIPHER)
library(Biostrings)
library(qdapRegex)
library_end_time <- Sys.time()
print(paste("Libraries loaded in:", difftime(library_end_time, library_start_time, units="mins"), "minutes"))

# Generate output files
output <- T

# List of bacteria runs
runs <- list.files(path.in)
print(paste("Found", length(runs), "runs to process:", paste(runs, collapse=", ")))

# Create folders
dir.create(path.rds, showWarnings=FALSE)
dir.create(path.out, showWarnings=FALSE)

# Settings
rc <- dada2:::rc
theme_set(theme_bw())

# Set DADA2 options to limit iterations and improve performance
print("Setting DADA2 options...")
dada2::setDadaOpt(MAX_CONSIST=10)
print(paste("MAX_CONSIST set to:", dada2::getDadaOpt("MAX_CONSIST")))

# ------------------------------------------------------------------------
# Loop over all bacteria runs
# ------------------------------------------------------------------------
for (run in runs) {
  run_start_time <- Sys.time()
  print(paste(run, "started at:", run_start_time))
  
  # ------------------------------------------------------------------------
  # DADA2 pipeline - PacBio
  # ------------------------------------------------------------------------
  
  # File names
  path_run <- paste(path.in, run, sep="")
  fns <- list.files(path_run, pattern="fastq", full.names=TRUE)
  print(paste("Found", length(fns), "fastq files in", run))
  assign(paste(run, "_fns", sep=""), fns)
  
  # Create directory for filtered files
  filter_dir <- file.path(path_run, "filtered")
  dir.create(filter_dir, showWarnings=FALSE, recursive=TRUE)
  
  # Sequence length before filtering
  print("Analyzing sequence lengths before filtering...")
  lens_start_time <- Sys.time()
  lens.fn <- lapply(fns, function(fn) nchar(getSequences(fn)))
  lens <- do.call(c, lens.fn)
  histo_unfiltered <- hist(lens, 100, main = paste(run, "unfiltered"))
  assign(paste(run, "histo_unfiltered", sep="_"), histo_unfiltered)
  lens_end_time <- Sys.time()
  print(paste("Length analysis completed in:", difftime(lens_end_time, lens_start_time, units="mins"), "minutes"))
  
  # Names filter
  filts <- file.path(path_run, "filtered", basename(fns))
  assign(paste(run, "_filts", sep=""), filts)
  
  # Filter
  print("Filtering and trimming sequences...")
  filter_start_time <- Sys.time()
  track <- filterAndTrim(fns, filts, minQ=3, minLen=500, maxLen=1800, maxN=0, rm.phix=T, maxEE=2, multithread=T)
  filter_end_time <- Sys.time()
  print(paste("Filtering completed in:", difftime(filter_end_time, filter_start_time, units="mins"), "minutes"))
  print(paste("Filtered from", sum(track[,1]), "to", sum(track[,2]), "reads"))
  
  assign(paste(run, "track", sep="_"), track)
  saveRDS(track, paste(path.rds, run, "_out.RDS", sep=""))
  
  # Sequence length after filtering
  histo_filtered <- hist(lens, 100, main = paste(run, "filtered"))
  assign(paste(run, "histo_filtered", sep="_"), histo_filtered)
  
  # Dereplicate
  print("Dereplicating sequences...")
  derep_start_time <- Sys.time()
  drp <- derepFastq(filts, verbose=TRUE)
  derep_end_time <- Sys.time()
  print(paste("Dereplication completed in:", difftime(derep_end_time, derep_start_time, units="mins"), "minutes"))
  
  assign(paste(run, "drp", sep="_"), drp)
  saveRDS(drp, paste(path.rds, run, "_derep.RDS", sep=""))
  
  # Error estimation for PacBio
  print("Learning error rates...")
  error_start_time <- Sys.time()
  err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=16, multithread=TRUE, nbases=1e8)
  error_end_time <- Sys.time()
  print(paste("Error learning completed in:", difftime(error_end_time, error_start_time, units="hours"), "hours"))
  
  assign(paste(run, "err", sep="_"), err)
  saveRDS(err, paste(path.rds, run, "_err.RDS", sep=""))
  
  # Denoise with inference algorithm
  print("Denoising sequences...")
  denoise_start_time <- Sys.time()
  dd <- dada(drp, err=err, BAND_SIZE=16, multithread=TRUE, selfConsist=TRUE)
  denoise_end_time <- Sys.time()
  print(paste("Denoising completed in:", difftime(denoise_end_time, denoise_start_time, units="hours"), "hours"))
  
  assign(paste(run, "dd", sep="_"), dd)
  saveRDS(dd, paste(path.rds, run, "_denoise.RDS", sep=""))
  
  run_end_time <- Sys.time()
  print(paste(run, "completed in:", difftime(run_end_time, run_start_time, units="hours"), "hours"))
}

# ------------------------------------------------------------------------
# Combine all runs
# ------------------------------------------------------------------------
print("Combining results from all runs...")
combine_start_time <- Sys.time()

# Combine and save
fns <- vector()
for(i in runs){
  fns <- c(fns, get(paste(i, "_fns", sep="")))
}
saveRDS(fns, paste(path.rds, "fns.RDS", sep=""))

filts <- vector()
for(i in runs){
  filts <- c(filts, get(paste(i, "_filts", sep="")))
}
saveRDS(filts, paste(path.rds, "filts.RDS", sep=""))

track <- get(paste(runs[1], "_track", sep=""))
if(length(runs)>1){
  for(i in runs[2:length(runs)]){
    track <- rbind(track, get(paste(i, "_track", sep="")))
  }
}
saveRDS(track, paste(path.rds, "track.RDS", sep=""))

drp <- vector()
for(i in runs){
  drp <- c(drp, get(paste(i, "_drp", sep="")))
}
saveRDS(drp, paste(path.rds, "drp.RDS", sep=""))

dd <- vector()
for(i in runs){
  dd <- c(dd, get(paste(i, "_dd", sep="")))
}
saveRDS(dd, paste(path.rds, "dd.RDS", sep=""))

combine_end_time <- Sys.time()
print(paste("Combination completed in:", difftime(combine_end_time, combine_start_time, units="mins"), "minutes"))

# Names from all bacterial runs
sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 2)

# Sequence table
print("Creating sequence table...")
seqtab_start_time <- Sys.time()
seqtab <- makeSequenceTable(dd)
print(paste("Sequence table dimensions:", paste(dim(seqtab), collapse=" x ")))
seqtab_end_time <- Sys.time()
print(paste("Sequence table creation completed in:", difftime(seqtab_end_time, seqtab_start_time, units="mins"), "minutes"))

# Remove chimeras
print("Removing chimeras...")
chimera_start_time <- Sys.time()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print(paste("Sequence table dimensions after chimera removal:", paste(dim(seqtab.nochim), collapse=" x ")))
print(paste("Percentage of sequences retained after chimera removal:", round(sum(seqtab.nochim)/sum(seqtab)*100, 2), "%"))
chimera_end_time <- Sys.time()
print(paste("Chimera removal completed in:", difftime(chimera_end_time, chimera_start_time, units="mins"), "minutes"))

# Tracking
getN <- function(x) sum(getUniques(x))
track <- cbind(track, sapply(dd, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nochim")
saveRDS(track, paste(path.rds, "track.RDS", sep=""))

# Taxa
print("Assigning taxonomy...")
taxa_start_time <- Sys.time()
taxa <- assignTaxonomy(seqtab.nochim, taxa_database, multithread=TRUE, tryRC=TRUE)
taxa_end_time <- Sys.time()
print(paste("Taxonomy assignment completed in:", difftime(taxa_end_time, taxa_start_time, units="mins"), "minutes"))

# Change headers from the long sequence to simple ASV numbers
seq <- colnames(seqtab.nochim)
ASV <- c()
for (i in 1:length(seq)) {
  ASV[i] <- paste("ASV", i, sep="")
}
ASV_seq <- data.frame(ASV=ASV, seq=seq)
ASV_seq[,1] <- paste(">", ASV_seq[,1], sep="")
saveRDS(ASV_seq, paste(path.rds, "ASV_seq.RDS", sep=""))

# Replace sequence headers by ASV numbers in seqtab.nochim
colnames(seqtab.nochim) <- ASV
saveRDS(seqtab.nochim, paste(path.rds, "seqtab.nochim.RDS", sep=""))

# Replace the sequence headers by the correct ASV number in taxa
for (i in 1:nrow(taxa)){
  for (j in 1:nrow(ASV_seq)){
    if(rownames(taxa)[i] == ASV_seq[j,2]){
      rownames(taxa)[i] <- ASV_seq[j,1]
    }
  }
}
rownames(taxa) <- gsub(">", "", rownames(taxa))

# "UTAX db" has different FASTA headers than "general release db"
# Taxonomy df needs to be done manually
# Empty df
temp <- data.frame(matrix(NA, nrow(taxa), 7))
rownames(temp) <- rownames(taxa)
colnames(temp) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

# Extract taxa
for(tax in 1:nrow(taxa)){
  tax_string <- paste0(taxa[tax,2],",") # Extract string and add "," at the end
  tax_vector <- qdapRegex::ex_between(tax_string, ":", ",")[[1]] # Use regex to extract taxa
  temp[tax, 1:length(tax_vector)] <- tax_vector # Add taxa
}

# Remove db-number
temp[,7] <- gsub("_SH[0-9]*.08FU", "", temp[,7])

# Overwrite taxa df
taxa <- temp; rm(temp)
saveRDS(taxa, paste(path.rds, "taxa.RDS", sep=""))

# ------------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------------
print("Starting sequence clustering...")
cluster_start_time <- Sys.time()

nproc <- 16

# Change ASV names to corresponding sequences
seqtab.nochim_ASV <- seqtab.nochim
colnames(seqtab.nochim_ASV) <- ASV_seq$seq

# Align sequences and compute distance matrix
print("Aligning sequences...")
align_start_time <- Sys.time()
ASV_seqs <- Biostrings::DNAStringSet(ASV_seq$seq) # Convert to DNA string set
ASV_seqs_aln <- DECIPHER::AlignSeqs(ASV_seqs, processors = nproc)
align_end_time <- Sys.time()
print(paste("Sequence alignment completed in:", difftime(align_end_time, align_start_time, units="mins"), "minutes"))

print("Computing distance matrix...")
dist_start_time <- Sys.time()
ASV_dist <- DECIPHER::DistanceMatrix(ASV_seqs_aln, processors = nproc)
dist_end_time <- Sys.time()
print(paste("Distance matrix computation completed in:", difftime(dist_end_time, dist_start_time, units="mins"), "minutes"))

# OTUs
print("Clustering into OTUs...")
otu_start_time <- Sys.time()
ASV97_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.03) # cutoff = 0.03 -> 97% OTU 
ASV98_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.02) # cutoff = 0.02 -> 98% OTU 
ASV99_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.01) # cutoff = 0.01 -> 99% OTU
otu_end_time <- Sys.time()
print(paste("OTU clustering completed in:", difftime(otu_end_time, otu_start_time, units="mins"), "minutes"))

# Loop over the 3 thresholds
for (n in c("97", "98", "99")) {
  threshold_start_time <- Sys.time()
  print(paste("Processing", n, "% OTUs..."))
  
  # Get the data
  ASV_cluster <- get(paste0("ASV", n, "_cluster"))
  ASV_cluster$ASV <- gsub(">","",ASV_seq$ASV) # Add ASV name
  ASV_cluster$ASV_seq <- ASV_seq$seq # Add sequences
  
  # Number of sequences per ASV
  print("Calculating sequence abundance per ASV...")
  ASV_cluster$ASV_abu <- NA
  for (seq in ASV_cluster$ASV_seq) {
    abu <- sum(seqtab.nochim_ASV[,colnames(seqtab.nochim_ASV) == seq])
    ASV_cluster$ASV_abu[ASV_cluster$ASV_seq == seq] <- abu
  }
  
  # Most abundant sequences for each OTU
  print("Finding most abundant sequences for each OTU...")
  ASV_cluster$OTU_seq <- NA
  ASV_cluster$OTU_abu <- NA
  
  for (OTU in ASV_cluster$cluster) {
    ASV_cluster_subset <- ASV_cluster[ASV_cluster$cluster == OTU,] # All entries for current OTU
    n_seq <- sum(ASV_cluster_subset$ASV_abu)
    most_abu_seq <- ASV_cluster_subset$ASV_seq[ASV_cluster_subset$ASV_abu == max(ASV_cluster_subset$ASV_abu)] # OTU(s) with most seqs
    most_abu_seq <- most_abu_seq[1] # Take the first entry in case multiple OTUs have max number of seqs
    ASV_cluster$OTU_seq[ASV_cluster$cluster == OTU] <- most_abu_seq # Add seq to table
    ASV_cluster$OTU_abu[ASV_cluster$cluster == OTU] <- n_seq # Add num to table
  }
  
  # Sort by OTU abundance
  ASV_cluster_sort <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),] # Sort by OTU name
  
  # Rename OTU by abundance
  clusters_sort <- unique(ASV_cluster_sort$cluster)
  ASV_cluster$OTU <- NA
  for(i in seq(length(clusters_sort))){
    ASV_cluster$OTU[ASV_cluster$cluster == clusters_sort[i]] <- paste0("OTU", i) # Rename
  }
  
  ### FINAL OTU FILES
  
  # OTU ASV_seq
  OTU_ASV_seq_all <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),]
  OTU_ASV_seq_all <- data.frame(OTU = OTU_ASV_seq_all$OTU, seq = OTU_ASV_seq_all$OTU_seq)
  OTU_ASV_seq <- data.frame(OTU = unique(OTU_ASV_seq_all$OTU), seq = NA)
  for (OTU in OTU_ASV_seq$OTU) {
    seq <- OTU_ASV_seq_all$seq[OTU_ASV_seq_all$OTU == OTU][1] # Get the sequence
    OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] <- seq # Add to seqtab.nochima frame
  }
  
  # OTU seqtab.nochim
  OTU_seqtab.nochim <- seqtab.nochim_ASV %>% t %>% rowsum(ASV_cluster$cluster) %>% t # Cluster
  for(i in seq(length(clusters_sort))){
    colnames(OTU_seqtab.nochim)[colnames(OTU_seqtab.nochim) == clusters_sort[i]] <- paste0("OTU", i) # Rename
  }
  OTU_seqtab.nochim <- OTU_seqtab.nochim[,OTU_ASV_seq$OTU] # Sort
  
  # OTU taxa
  OTU_taxa <- data.frame(matrix(nrow = nrow(OTU_ASV_seq), ncol = 7))
  rownames(OTU_taxa) <- OTU_ASV_seq$OTU
  colnames(OTU_taxa) <- colnames(taxa)
  
  for(OTU in OTU_ASV_seq$OTU){
    seq <- OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] # OTU seq
    ASV <- ASV_cluster$ASV[ASV_cluster$ASV_seq == seq] # ASV belong to OTU seq
    
    taxa_ASV <- taxa[rownames(taxa) == ASV[1], ] # Taxa belong to the first ASV
    OTU_taxa[rownames(OTU_taxa) == OTU, ] <- taxa_ASV # Transfer to OTU taxa
  }
  
  # OTU_ASV
  OTU_ASV <- data.frame(OTU = ASV_cluster$OTU, ASV = ASV_cluster$ASV,
                       OTU_seq = ASV_cluster$OTU_seq, ASV_seq = ASV_cluster$ASV_seq,
                       OTU_abu = ASV_cluster$OTU_abu, ASV_abu = ASV_cluster$ASV_abu)
  
  OTU_ASV <- OTU_ASV[order(ASV_cluster$OTU_abu, decreasing = T),] # Sort
  
  # Store clustered data
  assign(paste0("OTU_ASV_seq", n), OTU_ASV_seq)
  assign(paste0("OTU_seqtab.nochim", n), OTU_seqtab.nochim)
  assign(paste0("OTU_taxa", n), OTU_taxa)
  assign(paste0("OTU_ASV", n), OTU_ASV)
  
  # Save RDS files
  dir.create(paste0(path.rds,"ASV",n), showWarnings=FALSE)
  saveRDS(OTU_ASV_seq, paste0(path.rds,"ASV",n, "/SEQ", n, ".RDS"))
  saveRDS(OTU_seqtab.nochim, paste0(path.rds,"ASV",n, "/DAT", n, ".RDS"))
  saveRDS(OTU_taxa, paste0(path.rds,"ASV",n, "/TAXA", n, ".RDS"))
  saveRDS(OTU_ASV, paste0(path.rds,"ASV",n, "/ABU", n, ".RDS"))
  
  threshold_end_time <- Sys.time()
  print(paste(n, "% OTU processing completed in:", difftime(threshold_end_time, threshold_start_time, units="mins"), "minutes"))
}

cluster_end_time <- Sys.time()
print(paste("Clustering completed in:", difftime(cluster_end_time, cluster_start_time, units="hours"), "hours"))

# ------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------
if(output){
  output_start_time <- Sys.time()
  print("Generating output files...")
  
  # QUALITY CHECK
  dir.create(paste0(path.out,"quality_check"), showWarnings=FALSE)
  
  # Number of reads
  write.xlsx(track, paste0(path.out, "quality_check/track.xlsx"))
  
  # Length histo
  pdf(paste0(path.out,"quality_check/reads_length.pdf"))
  for(i in runs){
    plot(get(paste(i, "histo_unfiltered", sep="_")), main = paste(i, "unfiltered"))
    plot(get(paste(i, "histo_filtered", sep="_")), main = paste(i, "filtered"))
  }
  dev.off()
  
  # Quality profiles unfiltered and untrimmed
  pdf(paste0(path.out,"quality_check/reads_quality_unfilt_untrim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(fns[i]))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()
  
  # Quality profiles filtered and trimmed
  pdf(paste0(path.out,"quality_check/reads_quality_filt_trim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(filts[i]))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()
  
  # Error rate learning
  runs_err <- c(paste0(runs, "_err"))
  
  pdf(paste0(path.out,"quality_check/error_rate.pdf"))
  for (i in 1:length(runs)) {
    plot(plotErrors(get(runs_err[i]), nominalQ=TRUE))
    grid.text(runs[i], hjust=-1.5, vjust = -27.5, rot = 90)
  }
  dev.off()
  
  # Rarefaction
  pdf(paste0(path.out,"quality_check/rarefaction.pdf"))
  for (run in runs) {
    sp <- gsub("microbials", "m", run) # Search pattern
    DAT <- seqtab.nochim[grep(sp, rownames(seqtab.nochim)),]
    rarecurve(DAT, step = 20, label=F, main = paste("Rarefaction", run), ylab="ASVs", xlab="Sequencing depth")
  }
  dev.off()
  
  print("Quality check output done")
  
  # ASV TABLES
  
  # 97
  dir.create(paste0(path.out,"ASV97"), showWarnings=FALSE)
  write.table(OTU_seqtab.nochim97, paste0(path.out, "ASV97/DAT97.tab"), sep="\t")
  write.table(OTU_taxa97, paste0(path.out, "ASV97/TAXA97.tab"), sep="\t")
  write.table(OTU_ASV_seq97, paste0(path.out, "ASV97/SEQ97.tab"), sep="\t", row.names = F)
  write.table(OTU_ASV97, paste0(path.out, "ASV97/ABU.tab"), sep="\t", row.names = F)
  
  # 98
  dir.create(paste0(path.out,"ASV98"), showWarnings=FALSE)
  write.table(OTU_seqtab.nochim98, paste0(path.out, "ASV98/DAT98.tab"), sep="\t")
  write.table(OTU_taxa98, paste0(path.out, "ASV98/TAXA98.tab"), sep="\t")
  write.table(OTU_ASV_seq98, paste0(path.out, "ASV98/SEQ98.tab"), sep="\t", row.names = F)
  write.table(OTU_ASV98, paste0(path.out, "ASV98/ABU.tab"), sep="\t", row.names = F)
  
  # 99
  dir.create(paste0(path.out,"ASV99"), showWarnings=FALSE)
  write.table(OTU_seqtab.nochim99, paste0(path.out, "ASV99/DAT99.tab"), sep="\t")
  write.table(OTU_taxa99, paste0(path.out, "ASV99/TAXA99.tab"), sep="\t")
  write.table(OTU_ASV_seq99, paste0(path.out, "ASV99/SEQ99.tab"), sep="\t", row.names = F)
  write.table(OTU_ASV99, paste0(path.out, "ASV99/ABU.tab"), sep="\t", row.names = F)
  
  # 100
  dir.create(paste0(path.out,"ASV100"), showWarnings=FALSE)
  write.table(seqtab.nochim, paste0(path.out, "ASV100/DAT100.tab"), sep="\t")
  write.table(taxa, paste0(path.out, "ASV100/TAXA100.tab"), sep="\t")
  write.table(ASV_seq, paste0(path.out, "ASV100/SEQ100.tab"), sep="\t")
  
  output_end_time <- Sys.time()
  print(paste("Output file generation completed in:", difftime(output_end_time, output_start_time, units="mins"), "minutes"))
}

script_end_time <- Sys.time()
print(paste("Total script execution time:", difftime(script_end_time, script_start_time, units="hours"), "hours"))
