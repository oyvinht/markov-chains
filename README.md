# markov-chains

A Clojure library for working with Markov models.

## Usage

`(let [hmm (init-hmm {:s1 0.2 :s2 0.8}
                                        {:s1 {:s1 0.5 :s2 0.5}
                                         :s2 {:s1 0.3 :s2 0.7}}
                                        {:s1 {:N 0.3 :E 0.7}
                                         :s2 {:N 0.8 :E 0.2}})]
                      (viterbi (baum-welch hmm [:N :N :N :N :N :E :E :N :N :N] 1 [:N :E])
                               [:N :E :N]))
{:path (:s2 :s1 :s2), :prob 0.05031601966795285}`
