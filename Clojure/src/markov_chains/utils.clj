(ns markov-chains.utils)

(defn sumfn [func seq]
  "Apply func to all elements in seq and sum all results."
  (reduce (fn [res b] (+ res (func b))) 0 seq))

(defn n [m]
  "Helper: Normalize map of probabilities."
  (let [s (apply + (vals  m))]
    (if (== s 0)
      m
      (into {} (for [[k v] m] {k (/ v s)})))))

(defn m [seq-of-pairs]
  "Helper: Convert seq of 2-elem vectors into a map."
  (into {} seq-of-pairs))
