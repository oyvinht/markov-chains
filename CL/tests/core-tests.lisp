(eval-when (:compile-toplevel)
  (asdf:oos 'asdf:load-op 'fiveam))

(in-package :markov-chains)

(defmacro vec (&rest items)
  `(apply #'vector (quote ,items)))

#+5am
(5am:def-suite hmm-tests)

(5am:in-suite hmm-tests)

#+5am
(defvar *example-state-sequences*
  (list (vec a a a a b c a b d a d b f)
        (vec a a f a b b c d a a c)
        (vec a a b)
        (vec c a c f)
        (vec a a f b d e a)
        (vec f c c a b e d a a f a c)))

#+5am
(defvar *example-state-outcomes*
  (mapcar
   (lambda (el) (apply #'cons el))
   '((a 0)(a 0)(a 0)(c 0)(c 2)(b 0)(b 0)(e 0)(b 0)(a 1)(b 2)(a 1)(c 4)(d 1)(a 1)
     (b 0)(d 6)(e 0)(a 15)(c 4)(c 2)(a 0)(b 8)(a 2)(e 5)(a 2)(d 4)(c 20)(e 9))))

#+5am
(5am:test test-transition-probs
  "Checking probability of transitioning to state b when in a (25%)"
  (let ((hmm (make-hmm *example-state-sequences* *example-state-outcomes*)))
    (5am:is (= (transition-prob hmm 'a 'b) 0.25))))

#+5am
(5am:test test-outcome-probs
  "Checking probability for getting outcome 2 from state b (1/6)"
  (let ((hmm (make-hmm *example-state-sequences* *example-state-outcomes*)))
    (5am:is (= (outcome-prob hmm 'b 2) 1/6))))

;; Test forward algorithm using example from Russel & Norvig
#+5am
(defun make-example-init-probs-rn ()
  (let ((i (make-hash-table :test #'equalp)))
    (setf (gethash :rain i) 0.5)
    (setf (gethash :no-rain i) 0.5)
    i))

#+5am
(defun make-example-transition-probs-rn ()
  (let ((tp (make-hash-table :test #'equalp)))
    (let ((sub (make-hash-table :test #'equalp)))
      (setf (gethash :rain sub) 0.7)
      (setf (gethash :no-rain sub) 0.3)
      (setf (gethash :rain tp) sub))
    (let ((sub (make-hash-table :test #'equalp)))
      (setf (gethash :rain sub) 0.3)
      (setf (gethash :no-rain sub) 0.7)
      (setf (gethash :no-rain tp) sub))
    tp))

#+5am
(defun make-example-outcome-probs-rn ()
  (let ((op (make-hash-table :test #'equalp)))
    (let ((sub (make-hash-table :test #'equalp)))
      (setf (gethash :umbrella sub) 0.9)
      (setf (gethash :no-umbrella sub) 0.1)
      (setf (gethash :rain op) sub))
    (let ((sub (make-hash-table :test #'equalp)))
      (setf (gethash :umbrella sub) 0.2)
      (setf (gethash :no-umbrella sub) 0.8)
      (setf (gethash :no-rain op) sub))
    op))

#+5am
(5am:test test-hmm-rn-fw
          "Check if forward algorithm calculates correct alphas"
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
            (5am:is (= (forward hmm
                                (vector :umbrella :umbrella :no-umbrella :no-umbrella :umbrella)
                                2
                                :rain)
                       0.022965))))

#+5am
(5am:test test-hmm-rn-bw
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
             hmm (vector :umbrella :umbrella :no-umbrella :no-umbrella :umbrella))))
