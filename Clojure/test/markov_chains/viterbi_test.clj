(ns markov-chains.viterbi-test
  (:require [markov-chains.core :as c]
            [markov-chains.viterbi :as sut]
            [clojure.test :refer :all]))

(deftest test-hmm-viterbi-wp
  (testing "Testing Viterbi algorithm on Wikipedia example."
    (let [h (c/init-hmm {"Healthy" 0.6 "Fever" 0.4}
                        {"Healthy" {"Healthy" 0.7 "Fever" 0.3}
                         "Fever" {"Healthy" 0.4 "Fever" 0.6}}
                        {"Healthy" {"normal" 0.5 "cold" 0.4 "dizzy" 0.1}
                         "Fever" {"normal" 0.1 "cold" 0.3 "dizzy" 0.6}})]
      
      (is (= {:path '("Healthy" "Healthy" "Fever"), :prob 0.01512}
             (sut/viterbi h ["normal" "cold" "dizzy"]))))))
