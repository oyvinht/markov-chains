(in-package :markov-chains)

(defparameter *hash-table-test* (function equalp))

(defclass hmm ()
  ((transition-probabilities
    :accessor transition-probabilities
    :documentation "Per state probability of moving to each next state"
    :initarg :transition-probabilities
    :initform nil)
   (outcome-probabilities
    :accessor outcome-probabilities
    :documentation "Per state probability of observing one of M"
    :initarg :outcome-probabilities)
   (vocabulary
    :accessor vocabulary
    :documentation "Vocabulary of observations/emissions"
    :initarg :vocabulary)
   (init-probabilities
    :accessor init-probabilities
    :documentation "Per state probability of being the starting point"
    :initarg :init-probabilities)))

(defmethod initialize-instance :after ((hmm hmm) &key)
  (when (not (slot-boundp hmm 'vocabulary))
    (setf (vocabulary hmm) (states hmm))))

(defun make-hmm (init-probs-plist transition-probs-plist outcome-probs-plist)
  "Convenience function to make HMM from plists stating probabilities:

(make-hmm '(:s1 0.2 :s2 0.8)
          '(:s1 (:s1 0.5 :s2 0.5) :s2 (:s1 0.3 :s2 0.7))
          '(:s1 (:E 0.7 :N 0.3) :s2 (:E 0.2 :N 0.8)))"
  (let ((i-probs (make-hash-table))
        (t-probs (make-hash-table))
        (o-probs (make-hash-table)))
    (loop for (state prob) on init-probs-plist by #'cddr
          do (setf (gethash state i-probs) prob))
    (loop for (from-state transitions) on transition-probs-plist by #'cddr
          do (setf (gethash from-state t-probs)
                   (loop with to-state-probs = (make-hash-table)
                         for (to-state prob) on transitions by #'cddr
                         do (setf (gethash to-state to-state-probs) prob)
                         finally (return to-state-probs))))
    (loop for (state outcomes) on outcome-probs-plist by #'cddr
          do (setf (gethash state o-probs)
                   (loop with state-outcomes = (make-hash-table)
                         for (outcome prob) on outcomes by #'cddr
                         do (setf (gethash outcome state-outcomes) prob)
                         finally (return state-outcomes))))
    (make-instance 'hmm
                   :init-probabilities i-probs
                   :transition-probabilities t-probs
                   :outcome-probabilities o-probs)))

;;------------------------------------------------------------------------------
;; Various helper functions
;;------------------------------------------------------------------------------
(defun print-hashlist (hashlist)
  (loop for hash-table across hashlist
        for num from 1
        do (format t "~d.  " num)
           (maphash (lambda (state val)
                        (format t "~a  ~,8f  " state val))
                    hash-table)
           (format t "~%")))

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

(defun states (hmm)
  (loop for k being the hash-keys of (transition-probabilities hmm) collect k))

(defun outcome-probs (state-outcomes)
  "Create a nested hash table with outcome probabilities i->j->prob"
  (counts->probs (outcome-counts state-outcomes)))

(defun init-hmm (state-sequences state-outcomes)
  "Make a Hidden Markov Model from sets of initial state sequences and outcomes"
  (let ((tp (transition-probs state-sequences))
        (op (outcome-probs state-outcomes)))
    (make-instance
     'hmm :transition-probabilities tp :outcome-probabilities op)))

(defun init-prob (hmm state)
  (gethash state (init-probabilities hmm) 0))

(defun transition-prob (hmm src-state dst-state)
  "The probability of transitioning from src-state to dst-state for hmm"
  (gethash dst-state (gethash src-state (transition-probabilities hmm)) 0))

(defun outcome-prob (hmm state outcome)
  "The probability of observing outcome form state for hmm"
  (gethash outcome (gethash state (outcome-probabilities hmm)) 0))

;;------------------------------------------------------------------------------
;; Classic HMM functions
;;------------------------------------------------------------------------------
(defun forwards (hmm observations)
  "Build a vector of hash tables with forward propabilities (alphas) per state"
  (let ((alphas (make-array (length observations)
                            :initial-contents
                            (loop for num to (1- (length observations))
                                  collect (make-hash-table)))))
    ;; Base case
    (loop for state in (states hmm)
          do (setf (gethash state (aref alphas 0))
                   (* (init-prob hmm state)
                      (outcome-prob hmm state (aref observations 0)))))
    ;; Inductive cases
    (loop for obs-num from 1 to (1- (length observations))
          do (loop for state in (states hmm)
                   do (setf (gethash state (aref alphas obs-num))
                            (loop for prev-state in (states hmm)
                                  sum (* (gethash prev-state (aref alphas (1- obs-num)))
                                         (transition-prob hmm prev-state state)
                                         (outcome-prob hmm state (aref observations obs-num)))))))
    alphas))

(defun forward (hmm observations time state)
  "The forward probability of being in state at time given hmm and observations"
  (gethash state (aref (forwards hmm observations) time)))

(defun backwards (hmm observations)
  "Build a vector of hash tables with backward probabilities (betas) per state"
  (let ((betas (make-array (+ (length observations) 0)
                           :initial-contents
                           (loop for num to (- (length observations) 1)
                                 collect (make-hash-table)))))
    ;; Base case (always 100% since we know we're already there)
    (loop for state in (states hmm)
          do (setf (gethash state (aref betas (1- (length betas)))) 1))
    ;; Inductive cases (look at next last beta etc backwards)
    (loop for obs-num from (- (length observations) 2) downto 0
          for next-beta = (aref betas (1+ obs-num))
          do (loop for state in (states hmm)
                   do (setf
                       (gethash state (aref betas obs-num))
                       (loop for next-state being the hash-keys in next-beta
                             sum (* (gethash next-state next-beta)
                                    (transition-prob hmm state next-state)
                                    (outcome-prob hmm next-state (aref observations (1+ obs-num))))))))
    betas))

(defun backward (hmm observations time state)
  "The backward probability of being in state at time given hmm and observations"
  (gethash state (aref (backwards hmm observations) time)))

(defun test-hmm-rn-bw ()
  "Check if forward algorithm calculates correct betas"
  (let* ((hmm (make-instance 'hmm
                             :init-probabilities
                             (let ((htbl (make-hash-table)))
                               (setf (gethash :rain htbl) 0.5)
                               (setf (gethash :no-rain htbl) 0.5)
                               htbl)
                             :transition-probabilities
                             (let ((htbl (make-hash-table))
                                   (rain-htbl (make-hash-table))
                                   (no-rain-htbl (make-hash-table)))
                               ;; Prob of :rain/:no-rain when in :rain
                               (setf (gethash :rain rain-htbl) 0.7)
                               (setf (gethash :no-rain rain-htbl) 0.3)
                               (setf (gethash :rain htbl) rain-htbl)
                               ;; Prob of :rain/:no-rain when in :no-rain
                               (setf (gethash :rain no-rain-htbl) 0.3)
                               (setf (gethash :no-rain no-rain-htbl) 0.7)
                               (setf (gethash :no-rain htbl) no-rain-htbl)
                               htbl)
                             :outcome-probabilities
                             (let ((htbl (make-hash-table))
                                   (rain-htbl (make-hash-table))
                                   (no-rain-htbl (make-hash-table)))
                               ;; Prob of :umbrella/:no-umbrella when :rain
                               (setf (gethash :umbrella rain-htbl) 0.9)
                               (setf (gethash :no-umbrella rain-htbl) 0.1)
                               (setf (gethash :rain htbl) rain-htbl)
                               ;; Prob of :umbrella/:no-umbrella when :no-rain
                               (setf (gethash :umbrella no-rain-htbl) 0.2)
                               (setf (gethash :no-umbrella no-rain-htbl) 0.8)
                               (setf (gethash :no-rain htbl) no-rain-htbl)
                               htbl))))
    (backwards
     hmm (vector :umbrella :umbrella :no-umbrella :umbrella :umbrella))))
