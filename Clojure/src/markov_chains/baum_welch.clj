(ns markov-chains.baum-welch
  (:use [markov-chains.core]
        [markov-chains.utils]))

(defn- ξ [hmm observations t i j]
  (str "Expected state transition count:"
       "Probability of being in state 'i' at time 't'"
       "then in state 'j' at time 't+1'.")
  (when (not (get-in @(:cache hmm) [:ξ [observations t i j]]))
    (let [ksi (/ (* (forward hmm observations t i)
                    (transition-prob hmm i j)
                    (outcome-prob hmm j (nth observations (inc t)))
                    (backward hmm observations (inc t) j))
                 (reduce
                  (fn [res i]
                    (+ res
                       (reduce
                        (fn [res j]
                          (+ res
                             (* (forward hmm observations t i)
                                (transition-prob hmm i j)
                                (outcome-prob hmm j (nth observations (inc t)))
                                (backward hmm observations (inc t) j))))
                        0
                        (states hmm))))
                  0
                  (states hmm)))]
      (swap! (:cache hmm)
             (fn [cache] (assoc-in cache [:ξ [observations t i j]] ksi)))))
  (get-in @(:cache hmm) [:ξ [observations t i j]]))

(defn- γ [hmm observations t i]
  "Expected state occupancy."
  (when (not (get-in @(:cache hmm) [:γ [observations t i]]))
    (let [gamma (/ (* (forward hmm observations t i)
                      (backward hmm observations t i))
                   (reduce (fn [res s]
                             (+ res
                                (* (forward hmm observations t s)
                                   (backward hmm observations t s))))
                           0
                           (states hmm)))]
      (swap! (:cache hmm)
             (fn [cache] (assoc-in cache [:γ [observations t i]] gamma)))))
  (get-in @(:cache hmm) [:γ [observations t i]]))

(defn baum-welch [hmm observation-seqs num-iterations]
  "Return an improved HMM by training on a seq of observation-seqs."
  (letfn
      [(a [i j]
         ;; Excpected num transitions divided by total num transitions from i
         ;; for all observed sequences
         ;; NOTE: Parrallelizing heavy nominator
         (/ (apply + (map (fn [obs-seq]
                            (reduce (fn [res t] (+ res (ξ hmm obs-seq t i j)))
                                    0 (range (dec (count obs-seq)))))
                          observation-seqs))
            (reduce
             (fn [res obs-seq] (+ res(reduce
                                      (fn [res t] (+ res (γ hmm obs-seq t i)))
                                      0 (range (dec (count obs-seq))))))
             0 observation-seqs)))
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
           (reduce (fn [res s]
                     (conj res
                           {s (/ (sumfn (fn [obs-seq] (γ hmm obs-seq 0 s))
                                        observation-seqs)
                                 (count observation-seqs))}))
                   {}
                   (states hmm))
           ;; New state transition probabilities:
           (reduce (fn [res i]
                     (conj res
                           {i (n (reduce (fn [res j] (conj res {j (a i j)}))
                                         {}
                                         (states hmm)))}))
                   {}
                   (states hmm))
           ;; New state outcome probabilities
           (reduce
            (fn [res j] (conj res {j (n (reduce
                                         (fn [res o] (conj res {o (b j o)}))
                                         {} (:M hmm)))}))
            {} (states hmm)))]
      (if (<= num-iterations 1)
        new-hmm
        (recur new-hmm observation-seqs (dec num-iterations))))))
