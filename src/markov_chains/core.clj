(ns markov-chains.core)

;;;;---------------------------------------------------------------------------
;;;; General helper code.
;;;;---------------------------------------------------------------------------
(defn sum [fn seq]
  "Map fn over seq and sum all results."
  (apply + (map fn seq)))

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
  :π ; Per state probability of being the sequence starting point
  )

(defn make-hmm [state-sequences state-outcomes]
  "Create a HMM from a seq of state sequences and a seq of state outcomes."
  (let [tp (transition-probs state-sequences)
        op (outcome-probs state-outcomes)]
    (struct-map hmm :A (dissoc tp nil) :B op :π (get tp nil))))

(defn init-hmm [init-probs transition-probs outcome-probs]
  (struct-map hmm :A transition-probs :B outcome-probs :π init-probs))

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
  "Get probabilities of alle states at start time."
  (:π hmm))
       
(defn- n [m]
  "Helper: Normalize map of probabilities."
  (let [s (apply + (vals  m))]
    (if (== s 0)
      m
      (into {} (for [[k v] m] {k (/ v s)})))))

(defn- m [seq-of-pairs]
  "Helper: Convert seq of 2-elem vectors into a map."
  (into {} seq-of-pairs))

(defn forward [hmm observations]
  "Build an ordered list of maps with forward probabilities (alphas) per state."
  (loop [α [(m (map (fn [s] ; Base case
                      {s (* (init-prob hmm s)
                            (outcome-prob hmm s (first observations)))})
                    (states hmm)))]
         obs (rest observations)]
    (if (empty? obs)
      α
      (recur (concat α
                     [(m
                       (for [s_t+1 (states hmm)]
                         {s_t+1 (apply ; Inductive case
                                 + (map
                                    (fn [α_s]
                                      (* (val α_s)
                                         (transition-prob hmm (key α_s) s_t+1)
                                         (outcome-prob hmm s_t+1 (first obs))))
                                    (nth α (dec (count α)))))}))])
             (rest obs)))))

(defn sequence-prob [hmm observations]
  "The probability of seing a the sequence of observations given hmm: P(O|λ)."
  (let [fwds (forward hmm observations)]
    (sum (fn [s] (get (last fwds) s)) (states hmm))))

(defn backward [hmm observations]
  "Build an ordered list of maps with backward probabilities (betas) per state."
  (loop [β (list (m (map (fn [s] {s 1}) (states hmm))))
         obs observations]
    (if (empty? obs)
      β
      (recur (concat
              [(m (for [s (states hmm)]
                    {s (apply + (map (fn [β_s_t+1]
                                       (* (val β_s_t+1)
                                          (transition-prob hmm s (key β_s_t+1))
                                          (outcome-prob hmm
                                                        (key β_s_t+1)
                                                        (last obs))))
                                     (first β)))}))]
              β)
             (butlast obs)))))

(defn forward-backward [hmm observations]
  (let [fwds (forward hmm observations)
        bwds (backward hmm observations)]
    (map (fn [fw bw] (m (for [[k v] fw] {k (* v (get bw k))})))
         (cons (init-probs hmm) fwds)
         bwds)))

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
(defn rand-transition-probs [states]
  "Create a nested hash map with random transition probs between states."
  (into {}
        (map (fn [i]
               {i (into {} (map (fn [j] {j (rand)}) states))})
             states)))

(defn rand-outcome-probs [states vocabulary]
  "Create a nested hash map with random outcome probs from states."
  (into {}
        (map (fn [i]
               {i (into {} (map (fn [j] {j (rand)}) vocabulary))})
             states)))

(defn random-hmm [num-hidden-states]
  (let [init-states (take num-hidden-states (repeatedly (fn [] (gensym "S-"))))
        init (into {} (for [s init-states] {s (/ 1 num-hidden-states)}))
        A (rand-transition-probs init-states)
        B (rand-outcome-probs init-states '[N E])]
    (init-hmm init A B)))

(defn baum-welch [hmm observations num-iterations]
  "Return an improved HMM by training on observations num-iteration times."
  (let [fwds (forward hmm observations)
        bwds (rest (backward hmm observations))]
    (letfn
        [(ξ [t i j]
           ;; Expected state transition count:
           ;; Probability of being in state 'i' at time 't'
           ;; then in state 'j' at time 't+1'.
           (/ (* (get (nth fwds t) i)
                 (transition-prob hmm i j)
                 (outcome-prob hmm j (nth observations (inc t)))
                 (get (nth bwds (inc t)) j))
              (sum (fn [i]
                     (sum (fn [j]
                            (* (get (nth fwds t) i)
                               (transition-prob hmm i j)
                               (outcome-prob hmm j (nth observations (inc t)))
                               (get (nth bwds (inc t)) j)))
                          (states hmm)))
                   (states hmm))))
         (γ [t i]
           ;; Expected state occupancy count.
           (/ (* (get (nth fwds t) i)
                 (get (nth bwds t) i))
              (sum (fn [s]
                     (* (get (nth fwds t) s)
                        (get (nth bwds t) s)))
                   (states hmm))))
         (a [i j]
           ;; Excpected num transitions divided by total num transitions from i.
           (/ (sum (fn [t] (ξ t i j)) (range (dec (count observations))))
              (sum (fn [t] (γ t i)) (range (dec (count observations))))))
         (b [j v_k]
           ;; Expected num times in j observing v_k divided by
           ;; by tot expected num of times in j.
           (/ (sum (fn [t] (if (= v_k (nth observations t)) (γ t j) 0))
                   (range (count observations)))
              (sum (fn [t] (γ t j)) (range (count observations)))))]
      (let [new-init (into {} (map (fn [s] {s (γ 0 s)}) (states hmm)))
            new-A (into {} (map (fn [i]
                                  {i (n (m (map (fn [j]
                                                  {j (a i j)})
                                                (states hmm))))})
                                (states hmm)))
            new-B  (into {} (map (fn [j] {j (n (m (map (fn [o]
                                                         {o (b j o)})
                                                       (set observations))))})
                                 (states hmm)))
            new-hmm (init-hmm new-init new-A new-B)]
        (comment
          (doall
           (map (fn [x]
                  (println
                   (γ x (nth (states hmm) x))
                   "="
                   (sum (fn [j] (ξ x (nth (states hmm) x) j)) (states hmm)) "?"
                   (= (sum (fn [j] (ξ x (nth (states hmm) x) j)) (states hmm))
                      (γ x (nth (states hmm) x)))))
                (range (count (states hmm))))))
        ;; (println "Forwards:")
        ;; (println fwds)
        ;; (println "New PI:" new-init)
        ;; (println "New A:")
        ;; ;(clojure.pprint/pprint new-A)
        ;; (println "New B:")
        ;; ;(clojure.pprint/pprint new-B)
        ;; (println "FWDS" (apply + (vals (last fwds))))
        (if (= num-iterations 0)
          new-hmm
          (baum-welch new-hmm observations (dec num-iterations)))
        ))))
  
(defn test-baum-welch []
  (let [hmm (init-hmm
             {:s1 0.2 :s2 0.8}
             {:s1 {:s1 0.5 :s2 0.5}
              :s2 {:s1 0.3 :s2 0.7}}
             {:s1 {:N 0.3 :E 0.7}
              :s2 {:N 0.8 :E 0.2}})]

    (baum-welch hmm [:N :N :N :N :N :E :E :N :N :N] 1)))
