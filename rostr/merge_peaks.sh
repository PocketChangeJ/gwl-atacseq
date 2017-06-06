#RS requires peaks
#RS provides annotation
#RS widenode


# Collect list of files based on sample names
FILE_LIST=""
for SUBSAMPLE in $SAMPLELIST
do
	FILE_LIST=$FILE_LIST" "$DIR_OUTPUT/$SUBSAMPLE/$SUBSAMPLE"_peaks.narrowPeak"
done


# Combine samples using cat
cat $FILE_LIST > $FILE_OUTPUT.narrowPeak_cat.txt


# Cut out columns
cut -f 1-3,5 \
	$FILE_OUTPUT.narrowPeak_cat.txt \
		> $FILE_OUTPUT.narrowPeak.bed


# Nasty sort to force specific order of chromosomes
for CHR in `seq 1 22` X Y MT; 
do 
	grep -P "^$CHR |^$CHR\t|^chr$CHR |^chr$CHR\t" $FILE_OUTPUT.narrowPeak.bed | sort -k2,2n; 
done > $FILE_OUTPUT.narrowPeak_sort.bed


# Merge peaks
$TOOL_BEDTOOLS \
	merge \
	-i $FILE_OUTPUT.narrowPeak_sort.bed \
	-c 4 \
	-o count,mean \
		> $FILE_OUTPUT.narrowPeak_merge.bed


# Annotate peaks using R
$TOOL_RSCRIPT \
	$SCRIPT_ANNOPEAKS \
	$FILE_TSS \
	$FILE_OUTPUT.narrowPeak_merge.bed \
	$FILE_OUTPUT.narrowPeak_annot.bed
