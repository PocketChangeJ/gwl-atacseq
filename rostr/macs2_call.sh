
#RS provides peaks

$TOOL_MACS2 \
	callpeak \
	-t $FILE_INPUT \
	-f BAM \
	-g hs \
	--nomodel \
	--nolambda \
	--name $SAMPLE \
	--outdir $DIR_OUTPUT/$SAMPLE
