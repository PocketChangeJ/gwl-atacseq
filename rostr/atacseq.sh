# Notes
# Needs DESeq2,gplots installed in R
# Appears to be written originally for reference with no 'chr' in chromosome names

# Move to this node dir instead
DIR_NODES=$DIR_BASE/pipelines/atacseq/

# Could specify anything you want to look for in this value, uses 'find -name'
TEXT_INPUT='*.bam'

# Just something to start with
PIPELINE='macs2_call merge_peaks cover_peaks sam_idxstats calc_rpkm deseq2'

# Fix something specific for our HPC using SGE
SGE_PE='threaded'

# Override default dry run approach with sge
#SCHEDULER=sge
SCHEDULER=dry

# Modulestack
module load bed/bedtools/2.25.0 #/2.25.0 was default during development
module load R/3.3.3 #/3.2.2 was default and did not accept DESeq2 installation
module load sambamcram/samtools/1.3 #/1.3 was default during development

# We could use set -a to enable exporting all variables from this point on but
# explicitly using export does make you realize what you are doing I guess

# Scriptshack
export DIR_SCRIPTS=$DIR_NODES/scripts/
export SCRIPT_RPKM=$DIR_SCRIPTS/rpkm.R
export SCRIPT_ANNOPEAKS=$DIR_SCRIPTS/annotate.R
export SCRIPT_DESEQ2=$DIR_SCRIPTS/deseq2.R

# Toolshack
export TOOL_MACS2=/gnu/profiles/per-program/macs2/bin/macs2 
export TOOL_BEDTOOLS=bedtools
export TOOL_RSCRIPT=Rscript
export TOOL_SAMTOOLS=samtools

# Fileshack
export FILE_TSS=/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_data/TSS_Refseq_hg19.txt

# Scheduling information
ARGUMENTS_macs2_call='wtime:01,00,00 memory:8G'
ARGUMENTS_cover_peaks='wtime:01,00,00 memory:24G' # Turns out bedtools coverage eats RAM
ARGUMENTS_calc_rpkm='wtime:01,00,00 memory:8G'
ARGUMENTS_deseq2='wtime:04,00,00 memory:16G'

# Enable Differential Expression using DEseq2 using a file that names the controls
if [ -f "$FILE_CONTROLS" ]
then
	export FILE_CONTROLS=`realpath $FILE_CONTROLS`
else
	echo variable not pointing to existing file: FILE_CONTROLS
	exit
fi