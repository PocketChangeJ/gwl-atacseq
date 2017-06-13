(define-module (atacseq)
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

(define-public r-atacseq-scripts
  (package
    (name "r-atacseq-scripts")
    (version "1.0")
    (source (origin
              (method url-fetch)
              ;; TODO:  The repository is still private, so we cannot download
              ;; the source code from a public location.
              (uri (string-append
                    "file:///hpc/cog_bioinf/cuppen/personal_data/"
                    "rjanssen2/sources/gwl-atacseq-master.tar.gz"))
              (sha256
               (base32
                "0m3gilwz37d929cdi6a2li696za50wqb94lnnvr3sj8cy99awi74"))))
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
           (chdir "gwl-atacseq-master/rostr/scripts")
           (install-file "annotate.R" script-dir)
           (install-file "deseq2.R" script-dir)
           (install-file "rpkm.R" script-dir)))))
    (native-inputs
     `(("gzip" ,gzip)
       ("tar" ,tar)))
    (home-page #f)
    (synopsis "Additional scripts for the ATACseq pipeline.")
    (description "Additional scripts for the ATACseq pipeline.")
    (license #f)))

(define-public (call-peaks-for-sample sample run-path)
  (process
    (name (string-append "atacseq-" sample "-call-peaks"))
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
        (if (zero?
             (apply +
              (map
               (lambda (file-pair)
                 (let ((sample-name (car file-pair))
                       (file-name (cdr file-pair))
                       (output-path (string-append (getcwd) "/peaks")))

                   (unless (access? output-path F_OK)
                     (mkdir output-path))

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
               files))))))
   (synopsis "Call peaks using MACS2")
   (description "This process calls peaks for every sample in 'data-inputs'
using MACS2.")))

(define-public (merge-peaks-for-samples samples run-path)
  (process
   (name (string-append "atacseq-merge-peaks"))
   (version "1.0")
   (package-inputs
    `(("grep" ,grep)
      ("r" ,r)
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
            (bedtools (string-append
                       #$@(assoc-ref package-inputs "bedtools") "/bin/bedtools"))
            (cat (string-append
                  #$@(assoc-ref package-inputs "coreutils") "/bin/cat"))
            (cut (string-append
                  #$@(assoc-ref package-inputs "coreutils") "/bin/cut"))
            (sort (string-append #$@(assoc-ref package-inputs "coreutils")
                                 "/bin/sort")))
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
                 ;; Annotate peaks using R
                 (zero? (system* rscript annotate-script
                                 #$(assoc-ref data-inputs "tss-refseq")
                                 "narrowPeak_merge.bed"
                                 "narrowPeak_annot.bed")))
            #t    ; When everything executed just fine.
            #f))) ; When something went wrong.
   (synopsis "Merge peaks obtained from 'call-peaks'")
   (description "This process merges the peaks obtained using 'call-peaks'.")))

;;
;; XXX: The following processes are untested/WIP.
;;

(define (peak-coverage-for-samples samples run-path)
  (process
   (name "atacseq-peak-coverage")
   (version "1.0")
   (package-inputs
    `(("bedtools" ,bedtools)))
   (data-inputs
    `(("samples" . ,samples)))
   (run-time (complexity
              (space (gigabytes 24))
              (time (hours 1))))
   ;; TODO: Not sure what should be passed in the '-b' parameter.
   (procedure
    #~(let ((bedtools (string-append #$@(assoc-ref package-inputs "bedtools")
                                     "/bin/bedtools")))
        (if (zero?
             (apply + (map (lambda (sample)
                             (system (string-append bedtools " coverage "
                                                    "-a narrowPeak_annot.bed "
                                                    "-b " sample " -sorted "
                                                    "> merge_peak_cov.bed")))
                           '#$@(assoc-ref data-inputs "samples"))))
            #t
            #f)))
   (synopsis "")
   (description "")))

(define (samtools-idxstats-for-samples samples run-path)
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
        (if (zero?
             (apply + (map (lambda (sample)
                             (system (string-append samtools " idxstats "
                                                    sample
                                                    "> " sample ".idxstats.txt")))
                           '#$@(assoc-ref data-inputs "samples"))))
            #t
            #f)))
   (synopsis "Retrieve statistics from a BAM index file")
   (description "This process retrieves statistics from a BAM index file.")))
