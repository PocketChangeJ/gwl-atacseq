;; Copyright (C) 2017  Roel Janssen <R.R.E.Janssen-10@umcutrecht.nl>

;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.

;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.

(define-module (atacseq-workflow)
  #:use-module (guix processes)
  #:use-module (guix workflows)
  #:use-module (guix gexp)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages statistics)
  #:use-module (umcu packages python))

;; Merge all the coverage files to a count table.

;; Calculate RPKMs for each peak for each sample

;; Normalize the raw counts (coverage) and calculate differential expression using DESeq2
;;     maybe separate autosomal and all. Or normalize females vs females / males vs males.

;; Sample clustering and PCA

;; Peak Quality Control (# peaks / sample, # DE peaks / sample)

;; Calculate and normalize coverage in 100 bp bins for visualization

;; ----------------------------------------------------------------------------
;; PROCESS DEFINITIONS
;; ----------------------------------------------------------------------------

;; Quality control of ATAC-seq data (# mapped reads / sample, # reads /
;; chromosome / sample, # Y/X-chromosomal reads (gender check)
;;
;; ATAC_QC_V1.R
(define-public quality-control
  (process
   (name "quality-control")
   (version "1.0")
   (package-inputs
    `(("r" ,r)
      ("r-ggplot2" ,r-ggplot2)
      ("r-reshape2" ,r-reshape2)
      ("r-gridextra" ,r-gridextra)
      ("samtools" ,samtools)))
   ;; TODO: Find proper run-time resource settings.
   (run-time (complexity
              (space 1)
              (time 1)))
   (procedure #f)
   (synopsis "Quality control reporting for ATAC-seq data")
   (description "This process generates a quality control report that
includes statistics on the number of mapped reads per sample, number
of reads per chromosome, per sample, and a it includes a gender check.")))

(define-public call-peaks
  (process
   (name "call-peaks")
   (version "1.0")
   (package-inputs `(("python-macs2" ,python-macs2)))
   ;; TODO: Find proper run-time resource settings.
   (data-inputs '("sample1.bam" "sample2.bam"))
   (run-time (complexity
              (space 1)
              (time 1)))
   (procedure
    ;; #~(let ((files #$data-inputs)
    ;;         (macs2 (string-append
    ;;                 #$@(assoc-ref package-inputs "python-macs2")
    ;;                 "/bin/macs2")))
    ;;     (system* macs2 "callpeak"
    ;;              ;; ChIP-seq treatment file.
    ;;              "-t" file

    ;;              ;; File format.  We could use 'AUTO' here instead for automatic
    ;;              ;; detection of the input file type.
    ;;              "-f" "BAM"

    ;;              ;; Effective genome size; 'hs' is a shortcut for human (2.7e9)
    ;;              "-g" "hs"

    ;;              ;; Whether or not to build the shifting model.  If True, MACS
    ;;              ;; will not build model. by default it means shifting
    ;;              ;; size = 100, try to set extsize to change it.
    ;;              ;; DEFAULT: False.
    ;;              "--nomodel"

    ;;              ;; If True, MACS will use fixed background lambda as local
    ;;              ;; lambda for every peak region.  Normally, MACS calculates
    ;;              ;; a dynamic local lambda to reflect the local bias due to
    ;;              ;; potential chromatin structure.
    ;;              "--nolambda"

    ;;              ;; Experiment name.
    ;;              "-n"
    ;;              ;; TODO: Get the sample name from somewhere.
    ;;              (string-append sample "_peaks.txt")
    ;;              "--outdir" (string-append (getcwd) "/peaks")))
    #f)
   (synopsis "Call peaks using MACS2")
   (description "This process calls peaks for every sample in
'data-inputs' using MACS2.")))

;; Merge the peaks from each sample to one merged peak file using Bedtools merge.
;; Peak annotation (distance to TSS).
(define-public merge-peaks
  (process
   (name "merge-peaks")
   (version "1.0")
   (package-inputs `(("bedtools" ,bedtools)))
   ;; TODO: Find proper run-time resource settings.
   (run-time (complexity
              (space 1)
              (time 1)))
   (procedure #f)
   (synopsis "Merge peaks into one peak file")
   (description "This process nerges the peaks obtained using 'call-peaks'
for each sample to one merged peak file using Bedtools's 'merge' subcommand.")))

;; Calculate the coverage of the merged peaks for each sample using Bedtools coverage.
(define-public peaks-coverage
  (process
   (name "peaks-coverage")
   (version "1.0")
   (package-inputs `(("bedtools" ,bedtools)))
   (run-time (complexity
              (space 1)
              (time 1)))
   (procedure #f)
   (synopsis "Calculate coverage of a merged peaks file")
   (description "This process calculates the coverage of the merged peaks
for each sample using Bedtools's 'coverage' subcommand.")))

;; ----------------------------------------------------------------------------
;; WORKFLOW DEFINITION
;; ----------------------------------------------------------------------------
(define-public gwl-atacseq
  (workflow
   (name "gwl-atacseq")
   (version "1.0")
   (processes
    `(,quality-control
      ,call-peaks
      ,merge-peaks
      ,peaks-coverage))
   (restrictions
    `((,merge-peaks ,call-peaks)
      (,peaks-coverage ,merge-peaks)))
   (synopsis "")
   (description "")))
