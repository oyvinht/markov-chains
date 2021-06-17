(ns markov-chains.baum-welch
  (:use [markov-chains.core]
        [markov-chains.utils]))

(def ^:private ξ
  (memoize
   (fn [hmm observations t i j]
     (str
      "Expected state transition count:"
      "Probability of being in state 'i' at time 't'"
      "then in state 'j' at time 't+1'.")
     (/ (* (forward hmm observations t i)
           (transition-prob hmm i j)
           (outcome-prob hmm j (nth observations (inc t)))
           (backward hmm observations (inc t) j))
        (sumfn (fn [i]
                 (sumfn (fn [j]
                          (* (forward hmm observations t i)
                             (transition-prob hmm i j)
                             (outcome-prob hmm j (nth observations (inc t)))
                             (backward hmm observations (inc t) j)))
                        (states hmm)))
               (states hmm))))))

(def ^:private γ
  (memoize
   (fn [hmm observations t i]
     "Expected state occupancy."
     (/ (* (forward hmm observations t i)
           (backward hmm observations t i))
        (sumfn (fn [s]
                 (* (forward hmm observations t s)
                    (backward hmm observations t s)))
               (states hmm))))))

(defn baum-welch [hmm observation-seqs num-iterations]
  "Return an improved HMM by training on a seq of observation-seqs."
  (letfn
      [(a [i j]
         ;; Excpected num transitions divided by total num transitions from i
         ;; for all observed sequences
         ;; NOTE: Parrallelizing heavy nominator
         (/ (apply + (pmap (fn [obs-seq] (sumfn (fn [t] (ξ hmm obs-seq t i j))
                                                (range (dec (count obs-seq)))))
                           observation-seqs))
            (sumfn (fn [obs-seq] (sumfn (fn [t] (γ hmm obs-seq t i))
                                        (range (dec (count obs-seq)))))
                   observation-seqs)))
       (b [j v_k]
         ;; Expected num times in j observing v_k divided by
         ;; by tot expected num of times in j.
         (/ (sumfn (fn [obs-seq] (sumfn (fn [t] (if (= v_k (nth obs-seq t))
                                                  (γ hmm obs-seq t j)
                                                  0))
                                        (range (count obs-seq))))
                   observation-seqs)
            (sumfn (fn [obs-seq] (sumfn (fn [t] (γ hmm obs-seq t j))
                                        (range (count obs-seq))))
                   observation-seqs)))]
    (let [new-hmm
          (init-hmm
           ;; New state start probabilities:
           (m (map (fn [s] {s (/ (sumfn (fn [obs-seq] (γ hmm obs-seq 0 s))
                                        observation-seqs)
                                 (count observation-seqs))})
                   (states hmm)))
           ;; New state transition probabilities:
           (m (map (fn [i] {i (n (m (map (fn [j] {j (a i j)})
                                         (states hmm))))})
                   (states hmm)))
           ;; New state outcome probabilities
           (m (map (fn [j] {j (n (m (map (fn [o] {o (b j o)})
                                         (:M hmm))))})
                   (states hmm))))]
      (if (<= num-iterations 1)
        new-hmm
        (recur new-hmm observation-seqs (dec num-iterations))))))
