#RS requires annotation
#RS provides coverage

$TOOL_BEDTOOLS \
	coverage \
	-a $DIR_OUTPUT/wide/wide.narrowPeak_annot.bed \
	-b $FILE_INPUT \
	-sorted \
		> $FILE_OUTPUT.merged_peak_cov.bed
