(define-module (atacseq)
  #:use-module (ice-9 rdelim)
  #:use-module (guix processes)
  #:use-module (guix workflows)
  #:use-module (guix gexp)
  #:use-module (guix download)
  #:use-module (guix build-system trivial)
  #:use-module (guix packages)
  #:use-module (gnu packages base)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages statistics)
  #:use-module (umcu packages python)
  #:use-module (umcu packages fastqc))

;; ----------------------------------------------------------------------------
;; PACKAGES
;; ----------------------------------------------------------------------------

(define-public r-atacseq-scripts
  (let ((commit "56111981292d9f8eb66ba68bdaefd9d39b0f7f1e"))
    (package
     (name "r-atacseq-scripts")
     (version "1.0")
     (source (origin
              (method url-fetch)
              (uri (string-append
                    "https://github.com/UMCUGenetics/gwl-atacseq/archive/"
                    commit ".tar.gz"))
              (sha256
               (base32
                "0dpp4b5p03dmk0yp2710d1mpmiflfgl0s62gzyfhsfn2g64n5bzr"))))
     (build-system trivial-build-system)
     (arguments
      `(#:modules ((guix build utils))
        #:builder
        (begin
          (use-modules (guix build utils))
          (let* ((out (assoc-ref %outputs "out"))
                 (script-dir (string-append out "/share/atacseq/scripts/"))
                 (tar  (string-append (assoc-ref %build-inputs "tar") "/bin/tar"))
                 (PATH (string-append (assoc-ref %build-inputs "gzip") "/bin")))
            (mkdir-p script-dir)
            (setenv "PATH" PATH)
            (system* tar "xvf" (assoc-ref %build-inputs "source"))
            (chdir (string-append "gwl-atacseq-" ,commit "/rostr/scripts"))
            (install-file "annotate.R" script-dir)
            (install-file "deseq2.R" script-dir)
            (install-file "rpkm.R" script-dir)))))
     (native-inputs
      `(("gzip" ,gzip)
        ("tar" ,tar)))
     (propagated-inputs
      `(("r-genomicranges" ,r-genomicranges)))
     (home-page #f)
     (synopsis "Additional scripts for the ATACseq pipeline.")
     (description "Additional scripts for the ATACseq pipeline.")
     (license #f))))


;; ----------------------------------------------------------------------------
;; PROCESSES AND PROCESS TEMPLATE FUNCTIONS
;; ----------------------------------------------------------------------------

(define-public (call-peaks-for-sample sample run-path)
  (process
    (name (string-append sample "-call-peaks"))
    (version "1.0")
    (package-inputs `(("python-macs2" ,python-macs2)))
    (data-inputs
     `((,sample . ,(string-append
                   run-path "/" sample "/mapping/"
                   sample "_dedup.bam"))))
    (run-time (complexity
              (space (gigabytes 2))
              (time (minutes 20))))
   (procedure
    #~(let ((files '#$data-inputs)
            (macs2 (string-append #$@(assoc-ref package-inputs "python-macs2")
                                  "/bin/macs2")))
        (zero?
          (apply +
            (map
             (lambda (file-pair)
               (let ((sample-name (car file-pair))
                     (file-name (cdr file-pair))
                     (output-path (string-append (getcwd) "/peaks")))

                 (unless (access? output-path F_OK)
                   (catch #t
                     (lambda _ (mkdir output-path))
                     (lambda (key . arguments) #t)))

                 (system*
                  macs2 "callpeak"
                  ;; ChIP-seq treatment file.
                  "-t" file-name

                  ;; File format.  We could use 'AUTO' here instead for
                  ;; automatic detection of the input file type.
                  "-f" "BAM"

                  ;; Effective genome size; 'hs' is a shortcut for
                  ;; human (2.7e9)
                  "-g" "hs"

                  ;; Whether or not to build the shifting model.  If
                  ;; True, MACS will not build model. by default it means
                  ;; shifting size = 100, try to set extsize to change it.
                  ;; DEFAULT: False.
                  "--nomodel"

                  ;; If True, MACS will use fixed background lambda as
                  ;; local lambda for every peak region.  Normally, MACS
                  ;; calculates a dynamic local lambda to reflect the
                  ;; local bias due to potential chromatin structure.
                  "--nolambda"

                  ;; Experiment name.
                  "--name" sample-name
                  "--outdir" output-path)))
               files)))))
   (synopsis "Call peaks using MACS2")
   (description "This process calls peaks for every sample in 'data-inputs'
using MACS2.")))

(define-public (merge-peaks-for-samples samples)
  (process
   (name (string-append "merge-peaks"))
   (version "1.0")
   (package-inputs
    `(("grep" ,grep)
      ("r" ,r-minimal)
      ("bedtools" ,bedtools)
      ("coreutils" ,coreutils)
      ("atacseq-scripts" ,r-atacseq-scripts)))
   (data-inputs
    `(("samples" . ,samples)
      ("tss-refseq" . ,(string-append "/hpc/cog_bioinf/cuppen/project_data/"
                                      "Complex_svs/Common_data/"
                                      "TSS_Refseq_hg19.txt"))))
   (run-time (complexity
              (space (gigabytes 2))
              (time (hours 2))))
   (procedure
    #~(let ((egrep (string-append
                   #$@(assoc-ref package-inputs "grep") "/bin/egrep"))
            (rscript (string-append
                      #$@(assoc-ref package-inputs "r") "/bin/Rscript"))
            (annotate-script (string-append
                              #$@(assoc-ref package-inputs "atacseq-scripts")
                              "/share/atacseq/scripts/annotate.R"))
            ;; TODO: This doesn't take the propagated-inputs into account.
            (r-libs-site (string-append
                          #$@(assoc-ref package-inputs "atacseq-scripts")
                          "/site-library"))
            (bedtools (string-append
                       #$@(assoc-ref package-inputs "bedtools") "/bin/bedtools"))
            (cat (string-append
                  #$@(assoc-ref package-inputs "coreutils") "/bin/cat"))
            (cut (string-append
                  #$@(assoc-ref package-inputs "coreutils") "/bin/cut"))
            (sort (string-append #$@(assoc-ref package-inputs "coreutils")
                                 "/bin/sort")))

        ;; Remove files that will be overwritten.
        (when (access? "narrowPeak_sort.bed" F_OK)
          (delete-file "narrowPeak_sort.bed"))

        ;; Combine the samples
        (if (and (zero?
                  (apply +
                   (map (lambda (sample)
                          (let ((file (string-append (getcwd) "/peaks/" sample
                                                     "_peaks.narrowPeak")))
                            (system
                             (string-append cat " " file
                                            " >> narrowPeak_cat.txt"))))
                        '#$(assoc-ref data-inputs "samples"))))
                 ;; Cut out the useful columns
                 (zero? (system (string-append cut " -f 1-3,5 narrowPeak_cat.txt "
                                               "> narrowPeak.bed")))
                 ;; Nasty sort to force specific order of chromosomes
                 (zero?
                  (apply +
                   (map
                    (lambda (chr)
                      (let ((chr-str (cond
                                      ((number? chr) (number->string chr))
                                      ((symbol? chr) (symbol->string chr))
                                      ((string? chr) chr))))
                        (system
                         (string-append
                          egrep " \\\"^" chr-str "[[:space:]]|^chr" chr-str
                          "[[:space:]]\\\" narrowPeak.bed | "
                          sort " -k2,2n >> narrowPeak_sort.bed"))))
                    '(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
                      X Y MT))))
                 ;; Merge peaks
                 (zero?
                  (system
                   (string-append
                    bedtools " merge -i narrowPeak_sort.bed -c 4 -o count,mean "
                    "> narrowPeak_merge.bed")))
                 ;; Set the R environment so that GenomicRanges can be found.
                 ;; TODO: We need to fix the propagation of inputs here.
                 ;; FIXME: Yes, this is definitely cheating!
                 ;;(setenv "R_LIBS_SITE" r-libs-site)
                 (setenv "R_LIBS_SITE" "/gnu/profiles/per-language/r/site-library")
                 ;; Annotate peaks using R
                 (zero? (system* rscript annotate-script
                                 #$(assoc-ref data-inputs "tss-refseq")
                                 "narrowPeak_merge.bed"
                                 "narrowPeak_annot.bed")))
            (begin ; When everything executed just fine.
              (for-each delete-file '("narrowPeak.bed"
                                      "narrowPeak_cat.txt"
                                      "narrowPeak_sort.bed"
                                      "narrowPeak_merge.bed"))
              #t)
            #f))) ; When something went wrong.
   (synopsis "Merge peaks obtained from 'call-peaks'")
   (description "This process merges the peaks obtained using 'call-peaks'.")))

(define-public (peak-coverage-for-samples samples run-path)
  (process
   (name "peak-coverage")
   (version "1.0")
   (package-inputs
    `(("bedtools" ,bedtools)))
   (data-inputs
    `(("samples" . ,samples)))
   (run-time (complexity
              (space (gigabytes 48))
              (time (hours 1))))
   (procedure
    #~(let ((bedtools (string-append #$@(assoc-ref package-inputs "bedtools")
                                     "/bin/bedtools")))
        (zero?
          (apply + (map (lambda (sample)
                          (system
                           (string-append
                            ;; Run 'bedtools coverage'
                            bedtools " coverage "
                            ;; With the peak annotation information.
                            "-a narrowPeak_annot.bed "
                            ;; On the sample's BAM file.
                            "-b " #$run-path "/" sample "/mapping/" sample "_dedup.bam "
                            ;; And save the output to a file.
                            "-sorted > " sample "-merge_peak_cov.bed")))
                        '#$(assoc-ref data-inputs "samples"))))))
   (synopsis "Compute the coverage for each sample")
   (description "This process uses 'bedtools' to  compute both the depth and
breadth of coverage of features in each sample on the features in the union
of all provided samples. ")))

(define-public (peak-coverage-for-sample sample run-path)
  (process
   (name (string-append sample "-peak-coverage"))
   (version "1.0")
   (package-inputs
    `(("bedtools" ,bedtools)))
   (data-inputs sample)
   (run-time (complexity
              (space (gigabytes 70)) ; Peak usage was 58G
              (time (hours 1))))
   (procedure
    #~(let ((bedtools (string-append #$@(assoc-ref package-inputs "bedtools")
                                     "/bin/bedtools")))
        (zero? (system
                (string-append
                 ;; Run 'bedtools coverage'
                 bedtools " coverage "
                 ;; With the peak annotation information.
                 "-a narrowPeak_annot.bed "
                 ;; On the sample's BAM file.
                 "-b " #$run-path "/" #$data-inputs "/mapping/" #$data-inputs "_dedup.bam "
                 ;; And save the output to a file.
                 "-sorted > " #$data-inputs "-merge_peak_cov.bed")))))
   (synopsis (string-append "Compute the coverage for " sample))
   (description (string-append
                 "This process uses 'bedtools' to  compute both the depth and
breadth of coverage of features in " sample " on the features in the union
of the total set of provided samples."))))

(define-public (samtools-idxstats-for-samples samples run-path)
  (process
   (name "samtools-idxstats")
   (version "1.0")
   (package-inputs
    `(("samtools" ,samtools)))
   (data-inputs
    `(("samples" . ,samples)))
   (run-time (complexity
              (space (gigabytes 1))
              (time (minutes 20))))
   (procedure
    #~(let ((samtools (string-append #$@(assoc-ref package-inputs "samtools")
                                     "/bin/samtools")))
        (zero?
         (apply +
          (map (lambda (sample)
                 (system
                  (string-append
                   samtools " idxstats "
                   #$run-path "/" sample "/mapping/" sample "_dedup.bam "
                   "> " sample ".idxstats.txt")))
               '#$(assoc-ref data-inputs "samples"))))))
   (synopsis "Retrieve statistics from a BAM index file")
   (description "This process retrieves statistics from a BAM index file.")))

(define-public (calculate-rpkm-for-samples samples)
  (process
   (name "calculate-rpkms")
   (version "1.0")
   (package-inputs
    `(("r" ,r-minimal)
      ("atacseq-scripts" ,r-atacseq-scripts)))
   (data-inputs `(("samples" . ,samples)))
   (run-time (complexity
              (space (gigabytes 2))
              (time (hours 1))))
   (procedure
    #~(begin
        (use-modules (ice-9 format))
        (let ((rscript (string-append
                      #$@(assoc-ref package-inputs "r") "/bin/Rscript"))
            (rpkm-script (string-append
                          #$@(assoc-ref package-inputs "atacseq-scripts")
                          "/share/atacseq/scripts/rpkm.R")))
        ;; Create a sample list file.
        (with-output-to-file "samplelist.txt"
          (lambda _
            (for-each (lambda (sample)
                        (format #t "~a~/~a-merge_peak_cov.bed~/~a.idxstats.txt~%"
                                sample sample sample))
                      '#$(assoc-ref data-inputs "samples"))))
        ;; Run the rpkm.R script for the sample list.
        (zero? (system* rscript rpkm-script
                        "narrowPeak_annot.bed" "samplelist.txt"
                        (string-append (getcwd) "/RPKM"))))))
   (synopsis "Calculate RPKMs for samples")
   (description "This process merges the raw coverage files of each sample to
one count table and normalizes the coverage in ATAC-seq peaks using RPKMs.")))

(define-public (differential-expression controls-file)
  (process
   (name "differential-expression")
   (version "1.0")
   (package-inputs
    `(("r" ,r-minimal)
      ("atacseq-scripts" ,r-atacseq-scripts)))
   (data-inputs controls-file)
   (run-time (complexity
              (space (gigabytes 16))
              (time (hours 4))))
   (procedure
    #~(let ((rscript (string-append
                      #$@(assoc-ref package-inputs "r") "/bin/Rscript"))
            (deseq2-script (string-append
                          #$@(assoc-ref package-inputs "atacseq-scripts")
                          "/share/atacseq/scripts/deseq2.R")))
        ;; FIXME: This is cheating.
        (setenv "R_LIBS_SITE" "/gnu/profiles/per-language/r/site-library")
        (zero? (system* rscript deseq2-script
                        "RPKM.narrowPeak_annot_comb.bed"
                        #$data-inputs
                        ;; FIXME: Adjust deseq2.R for proper output paths.
                        (string-append (getcwd) "/DE")))))
   (synopsis "Differential expression")
   (description "This process performs a differential expression analysis
using DESeq2.")))


;; ----------------------------------------------------------------------------
;; WORKFLOW
;; ----------------------------------------------------------------------------
;;
;; The following construct enables us to hide all Scheme code from the end-user
;; by conditionally exposing symbols to the outside world.
;;
;; The condition we use here, is the availability of "samples.txt" in the
;; user's current working directory.

(define (samples-reader filename)
  "Returns a list of sample names read from FILENAME."

  (define (line-reader port items)
    (let ((line (read-line port)))
      (if (or (eof-object? line)
              (string= "" line))
          items
          (let ((sample (string-take line (string-index line #\tab))))
            (line-reader port (cons sample items))))))

  (call-with-input-file filename
    (lambda (port)
      ;; Skip the first line because it's a header.
      (read-line port)
      (line-reader port '()))))

(let ((samples-file (string-append (getcwd) "/samples.txt")))
  (when (access? (string-append (getcwd) "/samples.txt") F_OK)
    (let ((samples (samples-reader samples-file))
          (run-path (getcwd)))

      ;; Define the 'call-peaks' processes.
      (for-each (lambda (sample)
                  (primitive-eval
                   `(define-public
                      ,(symbol-append (string->symbol sample) '-call-peaks)
                      (call-peaks-for-sample ,sample ,run-path))))
                samples)

      ;; Define the 'peak-coverage' processes.
      (for-each (lambda (sample)
                  (primitive-eval
                   `(define-public
                      ,(symbol-append (string->symbol sample) '-peak-coverage)
                      (peak-coverage-for-sample ,sample ,run-path))))
                samples)

      ;; Define the 'merge-peaks' process
      (primitive-eval
       `(define-public merge-peaks (merge-peaks-for-samples ',samples)))

      ;; Define the 'calculate-rpkms' process
      (primitive-eval
       `(define-public calculate-rpkms (calculate-rpkm-for-samples ',samples)))

      ;; Define the 'samtools-idxstats' process
      (primitive-eval
       `(define-public idxstats (samtools-idxstats-for-samples ',samples ,run-path)))

      ;; Define the 'differential-expression' process
      (primitive-eval
       `(define-public diff-exp (differential-expression "sample-overview.txt")))

      (let ((peak-coverage-processes
             (map (lambda (sample)
                    (primitive-eval
                     (symbol-append (string->symbol sample) '-peak-coverage)))
                  samples))
            (call-peaks-processes
             (map (lambda (sample)
                    (primitive-eval (symbol-append (string->symbol sample)
                                                   '-call-peaks)))
                  samples)))
      (primitive-eval
       `(define-public atacseq-workflow
          ,(workflow
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
            (description ""))))))))
