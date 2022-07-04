(defpackage #:markov-chains
  (:use #:common-lisp #:cl-user)
  (:export #:backward
           #:forward
           #:make-hmm
           #:viterbi))
