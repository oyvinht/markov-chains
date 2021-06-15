(ns markov-chains.core)

;;;;---------------------------------------------------------------------------
;;;; General helper code.
;;;;---------------------------------------------------------------------------
(defn- sumfn [func seq]
  "Apply func to all elements in seq and sum all results."
  (reduce (fn [res b] (+ res (func b))) 0 seq))

(defn- n [m]
  "Helper: Normalize map of probabilities."
  (let [s (apply + (vals  m))]
    (if (== s 0)
      m
      (into {} (for [[k v] m] {k (/ v s)})))))

(defn- m [seq-of-pairs]
  "Helper: Convert seq of 2-elem vectors into a map."
  (into {} seq-of-pairs))

;;;;---------------------------------------------------------------------------
;;;; Code for generating matrices from example data.
;;;;---------------------------------------------------------------------------
(defn- transition-counts [state-sequences]
  "Create nested hash map with transition counts: i->j->count."
  (reduce
   (fn [hash-map state-sequence]
     (loop [states state-sequence
            i nil ; Previous state
            j (first states) ; Current state
            transition-counts hash-map]
       (if (not j)
         transition-counts
         (recur
          (rest states)
          j
          (first (rest states))
          (assoc-in transition-counts [i j]
                    (inc (or (get-in transition-counts [i j]) 0)))))))
   {}
   state-sequences))

(defn- outcome-counts [state-outcomes]
  "Create a nested hash map with with outcome M from state N: N->M->count."
  (reduce (fn [counts occurence]
            (assoc-in counts occurence (inc (or (get-in counts occurence) 0))))
          {}
          state-outcomes))

(defn- counts->probs [count-hash-map]
  "Convert from count hash-map to probability hash-map."
  (into {} (for [[k cnts] count-hash-map
                 :let [cnts-sum (apply + (vals cnts))]]
             {k (into {} (for [[v cnt] cnts]
                           {v (/ cnt cnts-sum)}))})))

(defn- transition-probs [state-sequences]
  "Create a nested hash map with transition probabilities i->j->prob."
  (counts->probs (transition-counts state-sequences)))

(defn- outcome-probs [state-outcomes]
  "Create a nested hash map with outcome probabilities k->v->prob."
  (counts->probs (outcome-counts state-outcomes)))

;;;;---------------------------------------------------------------------------
;;;; Code for creating HMMs.
;;;;---------------------------------------------------------------------------
(defstruct hmm
  :A ; Per state probability of moving to each next state
  :B ; Per state probability of observing any of M
  :M ; Vocabulary for observations/emissions
  :π) ; Per state probability of being the sequence starting point

(defn make-hmm [state-sequences state-outcomes]
  "Create a HMM from a seq of state sequences and a seq of state outcomes."
  (let [tp (transition-probs state-sequences)
        op (outcome-probs state-outcomes)]
    (struct-map
     hmm :A (dissoc tp nil) :B op :M (keys (first (vals op))) :π (get tp nil))))

(defn init-hmm [init-probs transition-probs outcome-probs]
  (struct-map hmm :A transition-probs :B outcome-probs
              :M (keys (first (vals outcome-probs))) :π init-probs))

;;;;---------------------------------------------------------------------------
;;;; Code for working on HMMs
;;;;---------------------------------------------------------------------------
(defn states [hmm]
  "Return a list of states (state names) for an HMM."
  (keys (:A hmm)))

(defn transition-prob [hmm from-state to-state]
  "Get probability of transitioning directly from from-state to to-state."
  (get-in (:A hmm) [from-state to-state] 0))

(defn outcome-prob [hmm state outcome]
  "Get probability of observing outcome from the given state."
  (get-in (:B hmm) [state outcome] 0))

(defn init-prob [hmm state]
  "Get probability of being in state at start time."
  (get (:π hmm) state 0))

(defn init-probs [hmm]
  "Get probabilities of all states at start time."
  (:π hmm))

(def ^:private forwards
  "Build an ordered list of maps with forward probabilities (alphas) per state."
  (memoize
   (fn [hmm observations]
     (reduce
      (fn [alphas o]
        (concat
         alphas
         [(reduce (fn [res s_t+1]
                    (conj
                     res {s_t+1 (sumfn ; Inductive case
                                 (fn [alpha_s]
                                   (* (val alpha_s)
                                      (transition-prob hmm (key alpha_s) s_t+1)
                                      (outcome-prob hmm s_t+1 o)))
                                 (last alphas))}))
                  {}
                  (states hmm))]))
      [(reduce ; Base case
        (fn [res s]
          (conj res
                {s (* (init-prob hmm s)
                      (outcome-prob hmm s (first observations)))}))
        {} (states hmm))]
      (rest observations)))))

(defn- forward [hmm observations time state]
  "Get the forward probability of state at time given hmm and observations."
  (get (nth (forwards hmm observations) time) state))

(defn sequence-prob [hmm observations]
  "The probability of seing the sequence of observations given hmm: P(O|λ)."
  (let [fwds (forwards hmm observations)]
    (sumfn (fn [s] (get (last fwds) s)) (states hmm))))

(def ^:private backwards
  "Build an ordered list of maps with backward probabilities (betas) per state."
  (memoize
   (fn [hmm observations]
     (reduce
      (fn [betas o]
        (cons (reduce
               (fn [res s]
                 (conj res
                       {s (sumfn (fn [beta_s_t+1]
                                   (* (val beta_s_t+1)
                                      (transition-prob hmm s (key beta_s_t+1))
                                      (outcome-prob hmm (key beta_s_t+1) o)))
                                 (first betas))}))
               {}
               (states hmm))
              betas))
      [(apply conj (map (fn [s] {s 1}) (states hmm)))]
      (reverse observations)))))

(defn- backward [hmm observations time state]
  "Get the backward probability of state at time given hmm and observations."
  (get (nth (rest (backwards hmm observations)) time) state))
  
(defn- argmax [paths]
  (loop [paths paths best (first paths)]
    (if (empty? paths)
      best
      (recur (rest paths)
             (if (> (:prob (first paths)) (:prob best)) (first paths) best)))))

(defn viterbi [hmm observations]
  "Find sequence of states with highest likelihood of explaining observations."
  (loop [obs (rest observations)
         paths (for [s (states hmm)]
                 {:path (list s)
                  :prob (* (init-prob hmm s)
                           (outcome-prob hmm s (first observations)))})]
    (if (empty? obs)
      (argmax paths)
      (recur
       (rest obs)
       (map (fn [state]
              (let [best (argmax
                          (map (fn [path]
                                 {:path (concat (:path path) (list state))
                                  :prob (* (:prob path)
                                           (transition-prob
                                            hmm (last (:path path)) state))})
                               paths))]
                {:path (:path best)
                 :prob (* (:prob best) (outcome-prob hmm state (first obs)))}))
            (states hmm))))))

;;;;---------------------------------------------------------------------------
;;;; Code for learning HMM from observations.
;;;;---------------------------------------------------------------------------
(defn- rand-transition-probs [states]
  "Create a nested hash map with random transition probs between states."
  (m (map (fn [i] {i (m (map (fn [j] {j (rand)}) states))}) states)))

(defn- rand-outcome-probs [states vocabulary]
  "Create a nested hash map with random outcome probs from states."
  (m (map (fn [i] {i (m (map (fn [j] {j (rand)}) vocabulary))}) states)))

(defn random-hmm [num-hidden-states vocabulary]
  "Create an HMM with random init, transition and outcome probabilities."
  (let [init-states (take num-hidden-states (repeatedly (fn [] (gensym "S-"))))
        init (m (for [s init-states] {s (/ 1 num-hidden-states)}))
        A (rand-transition-probs init-states)
        B (rand-outcome-probs init-states vocabulary)]
    (init-hmm init A B)))

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
