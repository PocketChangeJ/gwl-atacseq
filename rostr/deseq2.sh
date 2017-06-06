#RS requires annot_comb
#RS provides diffexp
#RS widenode

cp $FILE_CONTROLS $FILE_OUTPUT.controls

$TOOL_RSCRIPT \
	$SCRIPT_DESEQ2 \
	$FILE_OUTPUT.narrowPeak_annot_comb.bed \
	$FILE_OUTPUT.controls \
	$FILE_OUTPUT