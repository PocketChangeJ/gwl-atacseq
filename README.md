ATACseq pipeline
================

This repository contains an ATACseq pipeline.

To run the pipeline, you will need [GNU Guix](https://gnu.org/software/guix)
and the [Guix Workflow Language extension](https://git.roelj.com/guix/gwl).

Running the pipeline works like this:
```
git clone https://github.com/UMCUgenetics/gwl-atacseq.git
export GUIX_WORKFLOW_PATH=`pwd`/gwl-atacseq
guix workflow --run=gwl-atacseq
```
