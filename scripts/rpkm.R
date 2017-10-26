# This script merges the raw coverage files of each sample to one count table ("merged_peaks.bed")
# and is used to normalize (RPKMs) the coverage in ATAC-seq peaks

args             <- commandArgs(trailingOnly=TRUE)
file_mergedpeaks <- args[1]
file_sampleinfo  <- args[2]
file_output      <- args[3]

# Load files to work with
merged_peaks            <- read.delim (file_mergedpeaks, header = FALSE,
                                       stringsAsFactors = FALSE)

samplelist              <- read.delim (file_sampleinfo, header = FALSE,
                                       stringsAsFactors = FALSE)

# Setup header
row.names(merged_peaks) <- merged_peaks[,8]
merged_peaks            <- merged_peaks[,-8]
names(merged_peaks)     <- c("chr", "start", "end", "hits", "mean_score",
                             "nearest_gene", "distance_TSS")

# Merge all coverages into one peak table
for(i in 1:nrow(samplelist))
{
    row          <- samplelist[i,]
    input_file   <- read.delim (row[,2], header=FALSE, stringsAsFactors=FALSE)
    merged_peaks <- cbind (merged_peaks, input_file[,9])

    names(merged_peaks)[length(names(merged_peaks))] <- row[,1]
}

# Write this table to disk for the flying spaghetti monster knows why
# Perhaps we are meant to split this script into 2 scripts?
# Note: this used to overwrite the file it loaded previously into merged_peaks,
# which struck me as kind of weird behaviour.
write.table (merged_peaks,
             file = paste(file_output, ".narrowPeak_annot_comb.bed", sep = ""),
             quote = FALSE,
             sep = "\t",
             row.names = TRUE)

# Put statistics per sample into a dataframe
mapped_reads = data.frame()
for (i in 1:nrow(samplelist))
{
    row                       <- samplelist[i,]
    idxstats                  <- read.delim(row[,3], header= F)
    total_mapped_reads_sample <- sum(idxstats[,3])
    autosomal_reads_sample    <- sum(idxstats[idxstats[,1] %in% 1:22,3])

    rowtwo <- data.frame (sample          = row[,1],
                          mapped_reads    = total_mapped_reads_sample,
                          autosomal_reads = autosomal_reads_sample)

    mapped_reads              <- rbind(mapped_reads, rowtwo)
}

# Actually determine RPKMs for the dataset
merged_peaks_rpkms            <- merged_peaks[,c(1:7)]
merged_peaks_rpkms$peak_width <- merged_peaks_rpkms$end - merged_peaks_rpkms$start

for (i in 1:nrow(samplelist))
{
    row    <- samplelist[i,]
    sample <- row[,1]
    rpms   <- merged_peaks[,sample] / (mapped_reads[mapped_reads$sample == sample, "mapped_reads"] / 1e6)
    rpkms  <- rpms / (merged_peaks_rpkms$peak_width / 1000)

    merged_peaks_rpkms[,sample] <- rpkms
}

# And write RPKM data to file
write.table (merged_peaks_rpkms,
             file      = paste(file_output, ".rpkms.bed", sep = ""),
             quote     = FALSE,
             sep       = "\t",
             row.names = TRUE)
