(define-module (example-run)
  #:use-module (atacseq))

(define run-path "/path/to/processed_data/output")

;; Create a call-peaks process per sample.
(define-public sample1-call-peaks (call-peaks-for-sample "sample1" run-path))
(define-public sample2-call-peaks (call-peaks-for-sample "sample2" run-path))
(define-public sample3-call-peaks (call-peaks-for-sample "sample3" run-path))

;; Merge the peeks of all samples in a single process.
(define-public merge-peaks (merge-peaks-for-samples
                            '("sample1" "sample2" "sample3")
                            run-path))
