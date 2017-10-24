### 13-04-2017
library(DESeq2)
library(ggplot2)
library(gplots)
library(reshape2)
library(RColorBrewer)
library(genefilter)

args = commandArgs(trailingOnly=TRUE)
file_mergedpeaks = args[1]
file_output = args[2]
qc = TRUE
plots = TRUE

options(scipen = 999)

# Raw counts, not RPKMS!
merged_peaks = read.delim(file_mergedpeaks, header = T, stringsAsFactors = F, row.names = 1, check.names = F)


### DESEQ2

# Seems that most DE peaks are actually on X-chromosome between males and females.
merged_peaks_autosomal = merged_peaks[merged_peaks$chr != "X" & merged_peaks$chr !="Y" & merged_peaks$chr != "MT",]
peak_counts = merged_peaks_autosomal[,c(8:ncol(merged_peaks_autosomal))]
peak_counts_average = rowMeans(peak_counts)

# Removing peaks with extreme coverage:

peak_counts = peak_counts[peak_counts_average<1000,]
print(paste(
    "# Removed ", length(peak_counts_average[peak_counts_average > 1000]),
    " of ", nrow(merged_peaks_autosomal),
    " peaks because peak coverage is more than 1000 counts", sep = ""))

print("### Starting differential expression analysis using DESeq2 ###")

subjects <- as.factor(colnames(peak_counts))
subjects <- gsub(subjects, pattern = "-ATAC-\\d+", replacement = "")
unique_subjects <- unique(subjects)

for (subject in unique_subjects){
  print(paste("## Calculating differential expression sample ", i, " of ", length(subjects), ": ", subject, sep = "" ))

  condition <- ifelse(subjects == subject, subject, "control")
  condition <- factor(condition, levels=c("control", subject))
  colData <- data.frame(condition = condition, row.names = names(peak_counts))

  dds <- DESeqDataSetFromMatrix(countData = peak_counts,
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)

  write.table(as.data.frame(res), file = paste(DE_output,date, subject, "_DE_peaks.bed", sep = ""), quote = F, sep = "\t", row.names = T)

  MA_folder <- paste(DE_output, "MA-plots/", sep = "")
  ifelse(!dir.exists(file.path(MA_folder)), dir.create(file.path(MA_folder), showWarnings = FALSE), FALSE)

  pdf(file = paste(MA_folder, output_name, "_", subject,  ".pdf",sep = ""), width = 5, height = 5, pointsize = 10)
  plotMA(res, main=paste("MA-plot_", subject, sep = ""), ylim=c(-4,4),colNonSig = rgb(211,211,211,150, maxColorValue = 255), cex = 0.5)
  dev.off()

  res2 <- res[!is.na(res$padj), ]
  padj <- res2[res2$padj < 0.1,]
  pval <- res2[res2$pvalue < 0.05,]
  log2_neg <- res2[res2$log2FoldChange < -1,]
  log2_pos <- res2[res2$log2FoldChange > 1,]

  i <- i+1
}
normalized_counts = counts(dds, normalized = TRUE)
write.table(normalized_counts, file = paste(file_output, ".DE_peaks_normalized.txt", sep = ""),quote = F, sep = "\t", row.names = T)

#Pasting from RNAseq pipeline now

# In order to test for differential expression, we operate on raw counts and use discrete distributions.
# However for other downstream analyses ( e.g.  for visualization or clustering ) it
# might be useful to work with transformed versions of the count data.
# The function rlog, stands for regularized log, transforming the original count data to the log2 scale
# by fitting a model with a term for each sample and a prior distribution on the coefficients which is
# estimated from the data.

rld <- rlog(dds)
rldMat <- assay(rld)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#heatmap of the sample-to-sample distances
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)
pdf(paste(file_output, ".DE_distances.pdf",sep=""))
#pdf(file = paste(file_output, ".DE_peaks.pdf", sep = ""))
heatmap.2(mat,Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()
#principal component plot of the samples
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#png(paste(outdir,"$runname\_PCAplot.png",sep=""))
pdf(paste(file_output, ".DE_PCAplot.pdf",sep=""))
ggplot(data, aes(PC1,PC2,color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()


# add gene symbols to results table
# # # ----- GENE ANNOTATIONS --------
# # convert SYMBOL map to table
# s = toTable($shortname.egSYMBOL)
# # convert ENSEMBL map to table
# t = toTable($shortname.egENSEMBL)
# # merge symbol map table to ensembl map table
# gene_annotations = merge(t, s, by.x = "gene_id", by.y = "gene_id")

# resOrdered <- res[order(res$padj),] #order by adjusted pvalue
# resOrdered\$id = rownames(resOrdered)
# res1 = as.data.frame(resOrdered)
# res2 = merge(res1, gene_annotations, by.x = "id", by.y = "ensembl_id", all.x=TRUE)
# res4 = res2[order(res2\$pval),]
# res5 = res4[which(!duplicated(res4\$id)),]
# tablename <- paste(outdir,"$runname\_DEanalysis_all.txt",sep="")
# write.table(res5, file=tablename, sep = "\\t", col.names = T, row.names = F, quote = F)
