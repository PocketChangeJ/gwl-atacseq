#RS requires coverage idxstats
#RS provides rpkm annot_comb
#RS widenode

# Create a list of sample names and corresponding file paths
for SUBSAMPLE in $SAMPLELIST; 
do 
	echo -e $SUBSAMPLE"\t"$DIR_OUTPUT/$SUBSAMPLE/$SUBSAMPLE".merged_peak_cov.bed\t"$DIR_OUTPUT/$SUBSAMPLE/$SUBSAMPLE".idxstats.txt"
done > $FILE_OUTPUT.samplelist

$TOOL_RSCRIPT $SCRIPT_RPKM \
	$FILE_OUTPUT.narrowPeak_annot.bed \
	$FILE_OUTPUT.samplelist \
	$FILE_OUTPUT