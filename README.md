ATACseq pipeline
================

This repository contains an ATACseq pipeline.

To run the pipeline, you will need [GNU Guix](https://gnu.org/software/guix)
and the [Guix Workflow Language extension](https://git.roelj.com/guix/gwl).

Pipeline processes
------------------

This pipeline can do the following computation for you:

* Call peaks using MACS2
* Merge peaks called by MACS2
* Calculate the coverage per sample
* Normalize the coverage in ATAC-seq peaks
* Differential expression analysis using DESeq2


Running the workflow
--------------------
Running the pipeline works like this:

```
# Get the pipeline.
git clone https://github.com/UMCUgenetics/gwl-atacseq.git
export GUIX_WORKFLOW_PATH=`pwd`/gwl-atacseq
```

Create a file called `samples.txt` that contains a line for each
sample that you want to analyze:

```
sample	name	folder	identifier	gender	type
sample1	sample1	/path/to/IAP/output	LIMS_ID	M	C
sample2	sample2	/path/to/IAP/output	LIMS_ID	F	S
sample3	sample3	/path/to/IAP/output	LIMS_ID	M	S
```

When all of the above is done, there's only one final step to take:

```
# Run the pipeline.
cd /directory/containing/samples.txt/
guixr workflow --run=atacseq
```

Or to run it on a `grid-engine` cluster, replace the last command with:
```
guixr workflow --run=atacseq -e grid-engine
```

Each step can be run separately:
```
guixr process --list-available
guixr process -r <name-of-the-process>
```

A plot of the workflow can be generated:
```
guixr workflow --graph=atacseq | dot -Tpdf -o atacseq-workflow.pdf
```

License
-------

The code in this repository is licensed under the GPLv3.
