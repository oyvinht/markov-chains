(in-package :markov-chains)

(defparameter *hash-table-test* (function equalp))

(defclass hmm ()
  ((A
    :accessor hmm-a
    :documentation "Per state probability of moving to each next state"
    :initarg :A)
   (B
    :accessor hmm-b
    :documentation "Per state probability of observing one of M"
    :initarg :B)
   (M
    :accessor hmm-m
    :documentation "Vocabulary of observations/emissions"
    :initarg :M)
   (P
    :accessor hmm-p
    :documentation "Per state probability of being the starting point"
    :initarg :P)))

(defun counts->probs (counts)
  "Make hash table of probabilities, i->j->prob, from hash table of counts"
  (let ((probs (make-hash-table :test #'equalp)))
    (maphash
     #'(lambda (i js) ; Iterate over every i in counts
         (let ((total (loop for val being the hash-values of js summing val))
               (sub (make-hash-table :test #'equalp)))
           (maphash #'(lambda (j count) ; Calc % of total for each j
                        (setf (gethash j sub) (/ count total)))
                    js)
           (setf (gethash i probs) sub)))
     counts)
    probs))

(defun outcome-counts (state-outcomes)
  "Create a nested hash table pointing to outcome M from state N: N->M->count"
  (let ((counts (make-hash-table :test *hash-table-test*)))
    (loop for (N . M) in state-outcomes
          do (let ((counts-for-N (gethash N counts (make-hash-table
                                                    :test *hash-table-test*))))
               (incf (gethash M counts-for-N 0))
               (setf (gethash N counts) counts-for-N)))
    counts))

(defun transition-counts (state-sequences)
  "Create a nested hash table with counts of transitions i->j->counts"
  (let ((counts (make-hash-table :test *hash-table-test*)))
    (dolist (seq state-sequences)
      (loop
        for src = dst and dst across seq
        do (let ((sub (gethash
                       src counts (make-hash-table :test *hash-table-test*))))
             (incf (gethash dst sub 0))
             (setf (gethash src counts) sub))))
    counts))

(defun transition-probs (state-sequences)
  "Create a nested hash table with transition probabilities i->j->prob"
  (counts->probs (transition-counts state-sequences)))

(defun outcome-probs (state-outcomes)
  "Create a nested hash table with outcome probabilities i->j->prob"
  (counts->probs (outcome-counts state-outcomes)))

(defun make-hmm (state-sequences state-outcomes)
  "Make a Hidden Markov Model from sets of initial state sequences and outcomes"
  (let ((tp (transition-probs state-sequences))
        (op (outcome-probs state-outcomes)))
    (make-instance 'hmm :A tp :B op)))

(defmethod transition-prob ((hmm hmm) src-state dst-state)
  "The probability of transitioning from src-state to dst-state for hmm"
  (gethash dst-state (gethash src-state (hmm-a hmm))))

(defmethod outcome-prob ((hmm hmm) state outcome)
  "The probability of observing outcome form state for hmm"
  (gethash outcome (gethash state (hmm-b hmm))))
