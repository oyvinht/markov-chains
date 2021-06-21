(ns markov-chains.baum-welch-test
  (:require [clj-async-profiler.core :as prof]
            [markov-chains.baum-welch :as sut]
            [markov-chains.core :as c]
            [clojure.test :refer :all]))

(deftest test-ksi-gamma
  (testing "Seeing if gamma approx. equals sum of ksi for random HMM."
    (let [h (c/random-hmm 4 [:A :B :C])
          obs [:A :A :B :C :A]]
      (is
       (every?
        true?
        (map (fn [state-num]
               (= (float
                   (#'sut/γ h obs state-num (nth (c/states h) state-num)))
                  (float
                   (apply
                    + (map
                       (fn [j] (#'sut/ξ h
                                        obs
                                        state-num
                                        (nth (c/states h) state-num)
                                        j))
                       (c/states h))))))
             (range (count (c/states h)))))))))

(deftest test-baum-welch-1
  (testing "Testing if Baum-Welch converges towards correct hi-prec. value."
    (let [hmm (sut/baum-welch
               (c/init-hmm {:s1 1/5 :s2 4/5}
                             {:s1 {:s1 1/2 :s2 1/2}
                              :s2 {:s1 3/10 :s2 7/10}}
                             {:s1 {:N 3/10 :E 7/10}
                              :s2 {:N 4/5 :E 1/5}})
               (list [:N :N :N :N :N :E :E :N :N :N]) 1)]
      (is
       (= 3364526423555802/15688595904862135
          (c/transition-prob hmm :s2 :s1))))))

(deftest test-baum-welch-2
  (testing "Testing if Baum-Welch converges towards correct value."
    (let [hmm (c/init-hmm {:s1 0.2 :s2 0.8}
                            {:s1 {:s1 0.5 :s2 0.5}
                             :s2 {:s1 0.3 :s2 0.7}}
                            {:s1 {:N 0.3 :E 0.7}
                             :s2 {:N 0.8 :E 0.2}})]
      (is
       (= (float 0.14285715)
          (float
           (c/transition-prob 
            (sut/baum-welch hmm [[:N :N :N :N :N :E :E :N :N :N]] 1000)
            :s2 :s1)))))))

(deftest test-baum-welch-3
  (testing "Performance is reasonable for multiple sequences."
    (let [hmm (c/init-hmm {:s1 0.2 :s2 0.8}
                            {:s1 {:s1 0.5 :s2 0.5}
                             :s2 {:s1 0.3 :s2 0.7}}
                            {:s1 {:N 0.3 :E 0.7}
                             :s2 {:N 0.8 :E 0.2}})
          seqs (take 20 (repeatedly (fn [] (shuffle [:N :N :N :N :N :E :E :N :N :N]))))]
      (prof/profile
       ;(criterium.core/report-result
       ;(criterium.core/quick-benchmark
       (sut/baum-welch hmm seqs 100 )
       ;{})
       ))))
