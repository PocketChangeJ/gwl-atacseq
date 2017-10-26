suppressMessages (library (GenomicRanges))

# Turn args into variables
args              <- commandArgs(trailingOnly=TRUE)
file_tss          <- args[1]
file_bed          <- args[2]
file_ann          <- args[3]

# Load annotation data
TSS               <- read.delim(file_tss)
TSS_G             <- GRanges(seqnames = TSS$chromosome_name,
                             IRanges(start = TSS$transcript_start,
                                     end   = TSS$transcript_start + 1),
                             ensembl  = TSS$ensembl_gene_id)

# Load peak data
peaks             <- read.delim(file_bed, header = F)
peaks_G           <- GRanges(seqnames = peaks$V1,
                             IRanges(start = peaks$V2,
                                     end   = peaks$V3),
                             score    = peaks$V4)

# Find nearest annotations
distances         <- distanceToNearest(peaks_G, TSS_G, ignore.strand=TRUE)
nearest_genes     <- TSS[subjectHits(distances),]

# Bind annotations to peaks
peaks_output      <- cbind (peaks, nearest_genes$ensembl_gene_id)
peaks_output      <- cbind (peaks_output, mcols (distances)$distance)
peaks_output$name <- paste ("peak", row.names (peaks_output), sep = "")

# Save annotations in bed format
write.table (peaks_output,
             file = file_ann,
             quote = FALSE,
             sep = "\t",
             row.names = FALSE,
             col.names = FALSE)
