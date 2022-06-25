(ns markov-chains.viterbi
  (:use markov-chains.core))

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
