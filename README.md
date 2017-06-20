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

Now is the time to create a configuration file.
This looks like the following:

```
;; CONFIGURE THE PARTS IN THIS SECTION
;; -------------------------------------------------------------------------
;; Save this file as my-configuration.scm
;; To use a different file name, adjust the module name in the
;; next line accordingly.
(define-module (my-configuration)
  #:use-module (atacseq)
  #:use-module (guix workflows))

;; Provide the path where the processed data from the IAP can be found.
(define run-path "/path/to/IAP/output")

;; Provide a list containing the samples that can be found in
;; the 'run-path' and that you want to analyze.
(define samples '("sample1" "sample2" "sample3"))

;; DO NOT EDIT ANYTHING BELOW THIS LINE.
;; -------------------------------------------------------------------------

;; Define a 'call-peaks' process for each sample.
(for-each (lambda (sample)
            (primitive-eval
             `(define-public
                ,(symbol-append (string->symbol sample) '-call-peaks)
                (call-peaks-for-sample ,sample run-path))))
          samples)

;; Do the same thing for 'peak-coverage'.
(for-each (lambda (sample)
            (primitive-eval
             `(define-public
                ,(symbol-append (string->symbol sample) '-peak-coverage)
                (peak-coverage-for-sample ,sample run-path))))
          samples)

;; These processes operate on all samples, so one for each is enough.
(define-public merge-peaks (merge-peaks-for-samples samples))
(define-public calculate-rpkms (calculate-rpkm-for-samples samples))
(define-public idxstats (samtools-idxstats-for-samples samples run-path))
(define-public diff-exp (differential-expression "sample-overview.txt"))

(define peak-coverage-processes
  (map (lambda (sample)
         (primitive-eval (symbol-append (string->symbol sample)
                                        '-peak-coverage)))
       samples))

(define call-peaks-processes
  (map (lambda (sample)
         (primitive-eval
          (symbol-append (string->symbol sample)
                         '-call-peaks)))
       samples))

(define-public atacseq-workflow
  (workflow
   (name "atacseq")
   (version "1.0")
   (processes
    (append `(,merge-peaks
              ,calculate-rpkms
              ,idxstats
              ,diff-exp)
            call-peaks-processes
            peak-coverage-processes))
   (restrictions
    (append
     ;; Before we can run 'merge-peaks', we need to 'call-peaks' on each sample.
     `((,merge-peaks . ,call-peaks-processes))
     ;; Each 'peak-coverage' process depends on the 'merge-peaks' process.
     (map (lambda (p) `(,p ,merge-peaks)) peak-coverage-processes)
     ;; Before we can 'calculate-rpkms', we need the output of all
     ;; 'peak-coverage' processes.
     `((,calculate-rpkms . ,(append peak-coverage-processes `(,idxstats))))
     `((,diff-exp ,calculate-rpkms))))
   (synopsis "")
   (description "")))

```

Save this file to the directory where the pipeline should place its output.
Now create a file called `sample-overview.txt` that contains the following
line for each sample that you want to analyze:
```
sample	name	folder	lims	sex
sample1	sample1	/path/to/IAP/output	LIMS_ID	M
sample2	sample2	/path/to/IAP/output	LIMS_ID	F
sample3	sample3	/path/to/IAP/output	LIMS_ID	M
```

When all of the above is done, there's only one final easy step to take:

```
# Run the pipeline.
cd /directory/containing/my-configuration.scm/
guixr workflow --run=atacseq
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
