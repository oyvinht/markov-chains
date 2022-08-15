(in-package :markov-chains)

(defun gamma (hmm observations timestep i)
  "Expected state occupancy"
  (/ (* (forward hmm observations timestep i)
        (backward hmm observations timestep i))
     (reduce (lambda (res s)
               (+ res
                  (* (forward hmm observations timestep s)
                     (backward hmm observations timestep s))))
             (states hmm)
             :initial-value 0)))

(defun ksi (hmm observations ts i j)
  "Expected state transition count:
Probability of being in state I at time TS then J at TS+1"
  (/ (* (forward hmm observations ts i)
        (transition-prob hmm i j)
        (outcome-prob hmm j (elt observations (1+ ts)))
        (backward hmm observations (1+ ts) j))
     (loop for i-st in (states hmm)
           sum (loop for j-st in (states hmm)
                     sum (* (forward hmm observations ts i-st)
                            (transition-prob hmm i-st j-st)
                            (outcome-prob hmm j-st (elt observations (1+ ts)))
                            (backward hmm observations (1+ ts) j-st))))))

(defun normalize (hash-table)
  "Return a new version of HASH-TABLE where the sum of values = 1"
  (let ((sum (loop for v being the hash-values of hash-table sum v)))
    (if (or (= sum 0) (= sum 1))
        hash-table
        (let ((new-hash-table (make-hash-table)))
          (loop for k being the hash-keys of hash-table
                do (setf (gethash k new-hash-table)
                         (/ (gethash k hash-table) sum))
                finally (return new-hash-table))))))

(defun baum-welch (hmm observation-sequences num-iterations)
  "Improve HMM by training on OBSERVATION-SEQUENCES"
  (declare (ignore num-iterations))
  (flet ((a (hmm i j)
           (/ (loop for obs-seq across observation-sequences
                    sum (loop
                          for ts from 0 to (- (length obs-seq) 2)
                          sum (ksi hmm obs-seq ts i j)))
              (loop for obs-seq across observation-sequences
                    for ts from 0
                    sum (gamma hmm obs-seq ts i))))
         (b (j v_k)
           (/ (loop for obs-seq across observation-sequences
                    sum (loop for ts from 0 to (- (length obs-seq) 1)
                              sum (if (eql v_k (elt obs-seq ts))
                                      (gamma hmm obs-seq ts j)
                                      0)))
              (loop for obs-seq across observation-sequences
                    sum (loop for ts from 0 to (- (length obs-seq) 1)
                              sum (gamma hmm obs-seq ts j))))))
    (psetf
     ;; New state start probabilities
     (init-probabilities hmm)
     (loop for state in (states hmm)
           with new-init-probs = (make-hash-table)
           do (setf (gethash state new-init-probs)
                    (/ (loop for obs-seq across observation-sequences
                             sum (gamma hmm obs-seq 0 state))
                       (length observation-sequences)))
           finally (return new-init-probs))
     ;; New state transition probabilities
     (transition-probabilities hmm)
     (loop for i in (states hmm)
           with i-probs = (make-hash-table)
           do (setf (gethash i i-probs)
                    (normalize
                     (loop for j in (states hmm)
                           with j-probs = (make-hash-table)
                           do (setf (gethash j j-probs)
                                    (a hmm i j))
                           finally (return j-probs))))
           finally (return i-probs))
     ;; New state outcome probabilities
     (outcome-probabilities hmm)
     (with-slots (outcome-probabilities) hmm
       (let ((new-outcome-probabilities (make-hash-table)))
         (loop for j in (states hmm)
               do (setf (gethash j new-outcome-probabilities)
                        (identity;normalize
                         (loop for o being the hash-keys of (gethash j outcome-probabilities)
                               with j-outcome-probs = (make-hash-table)
                               do (setf (gethash o j-outcome-probs)
                                        (b j o))
                               finally (return j-outcome-probs)))))
         new-outcome-probabilities)))
    hmm))
