(ns markov-chains.core-test
  (:require [clojure.test :refer :all]
            [markov-chains.core :as sut]))

;;;;---------------------------------------------------------------------------
;;;; Example state sequences and state outcomes.
;;;;---------------------------------------------------------------------------
(def example-state-sequences-1
  '[[a a a a b c a b d a d b f]
    [a a f a b b c d a a c]
    [a a b]
    [c a c f]
    [a a f b d e a]
    [f c c a b e d a a f a c]])

(def example-state-outcomes-1
  '[[a 0][a 0][a 0][c 0][c 2][b 0][b 0][e 0][b 0][a 1][b 2][a 1][c 4][d 1][a 1]
    [b 0][d 6][e 0][a 15][c 4][c 2][a 0][b 8][a 2][e 5][a 2][d 4][c 20][e 9]])

(deftest make-hmm-test
  (testing "Creating HMM from example data."
    (def hmm (sut/make-hmm example-state-sequences-1
                           example-state-outcomes-1))
    (is (not (nil? hmm)))))

(deftest transition-prob
  (testing "Checking probability of transitioning from state a to b (25%)."
    (let [hmm (sut/make-hmm example-state-sequences-1
                            example-state-outcomes-1)]
      (is (= (sut/transition-prob hmm 'a 'b) 1/4)))))

(deftest outcome-prob
  (testing "Checking probability of getting outcome 2 from state b (1/6)."
    (let [hmm (sut/make-hmm example-state-sequences-1
                            example-state-outcomes-1)]
      (is (= (sut/outcome-prob hmm 'b 2) 1/6)))))

;;;;---------------------------------------------------------------------------
;;;; Umbrella example from Russel & Norvig.
;;;;---------------------------------------------------------------------------
(def example-hmm-rn
  (sut/init-hmm {:rain 0.5 :no-rain 0.5}
                {:rain {:rain 0.7 :no-rain 0.3}
                 :no-rain {:rain 0.3 :no-rain 0.7}}
                {:rain {:umbrella 0.9 :no-umbrella 0.1}
                 :no-rain {:umbrella 0.2 :no-umbrella 0.8}}))

;; (deftest test-hmm-rn-fw
;;   (testing "Checking alphas of forward algorithm from Russel & Norvig."
;;     (is
;;      (= (list
;;          {:rain 0.81818181818181810, :no-rain 0.18181818181818182}
;;          {:rain 0.88335704125177800, :no-rain 0.11664295874822190}
;;          {:rain 0.19066793972352525, :no-rain 0.80933206027647480}
;;          {:rain 0.73079400458498200, :no-rain 0.26920599541501794}
;;          {:rain 0.86733888957548470, :no-rain 0.13266111042451528})
;;         (#'sut/forwards example-hmm-rn
;;                         [:umbrella :umbrella :no-umbrella :umbrella :umbrella])))))

;; (deftest test-hmm-rn-bw
;;   (testing "Testing backward algorithm on Russel & Norvig example."
;;     (is
;;      (= (list
;;          {:rain 0.6469355558301939, :no-rain 0.3530644441698061}
;;          {:rain 0.5923176018339928, :no-rain 0.4076823981660072}
;;          {:rain 0.37626717588941005, :no-rain 0.6237328241105898}
;;          {:rain 0.6533428165007112, :no-rain 0.34665718349928876}
;;          {:rain 0.6272727272727272, :no-rain 0.37272727272727274}
;;          {:rain 1, :no-rain 1})
;;         (#'sut/backwards
;;          example-hmm-rn
;;          [:umbrella :umbrella :no-umbrella :umbrella :umbrella])))))

;; (deftest test-hmm-rn-fw-bw
;;   (testing "Testing forward-backward on Russel & Norvig example."
;;     (is
;;      (= (list
;;          {:rain 0.64693555583019390 :no-rain 0.35306444416980610}
;;          {:rain 0.86733888957548470 :no-rain 0.13266111042451528}
;;          {:rain 0.82041905362367540 :no-rain 0.17958094637632463}
;;          {:rain 0.30748357600661774 :no-rain 0.69251642399338230}
;;          {:rain 0.82041905362367530 :no-rain 0.17958094637632466}
;;          {:rain 0.86733888957548470 :no-rain 0.13266111042451528})
;;         (sut/forward-backward
;;          example-hmm-rn
;;          [:umbrella :umbrella :no-umbrella :umbrella :umbrella])))))

;;;;---------------------------------------------------------------------------
;;;; Example from Wikipedia.
;;;;---------------------------------------------------------------------------
;; (deftest test-hmm-wp
;;   (testing "Testing forward-backward from Wikipedia example."
;;     (let [h (sut/init-hmm {"Healthy" 0.6 "Fever" 0.4}
;;                           {"Healthy" {"Healthy" 0.69 "Fever" 0.3 "E" 0.01}
;;                            "Fever" {"Healthy" 0.4 "Fever" 0.59 "E" 0.01}}
;;                           {"Healthy" {"normal" 0.5 "cold" 0.4 "dizzy" 0.1}
;;                            "Fever" {"normal" 0.1 "cold" 0.3 "dizzy" 0.6}})]
                          
;;       (is (= (list
;;               {"Healthy" 0.68308190362154340 "Fever" 0.31691809637845664}
;;               {"Healthy" 0.87701103755732590 "Fever" 0.12298896244267407}
;;               {"Healthy" 0.62322803095095390 "Fever" 0.37677196904904610}
;;               {"Healthy" 0.21095270484130565 "Fever" 0.78904729515869430})
;;              (sut/forward-backward h ["normal" "cold" "dizzy"]))))))

(deftest test-hmm-viterbi-wp
  (testing "Testing Viterbi algorithm on Wikipedia example."
    (let [h (sut/init-hmm {"Healthy" 0.6 "Fever" 0.4}
                          {"Healthy" {"Healthy" 0.7 "Fever" 0.3}
                           "Fever" {"Healthy" 0.4 "Fever" 0.6}}
                          {"Healthy" {"normal" 0.5 "cold" 0.4 "dizzy" 0.1}
                           "Fever" {"normal" 0.1 "cold" 0.3 "dizzy" 0.6}})]
                          
      (is (= {:path '("Healthy" "Healthy" "Fever"), :prob 0.01512}
             (sut/viterbi h ["normal" "cold" "dizzy"]))))))

;;;;---------------------------------------------------------------------------
;;;; Test Baum-Welch
;;;;---------------------------------------------------------------------------
(deftest test-ksi-gamma
  (testing "Seeing if gamma approx. equals sum of ksi for random HMM."
    (let [h (sut/random-hmm 4 [:A :B :C])
          obs [:A :A :B :C :A]]
      (is
       (every?
        true?
        (map (fn [state-num]
               (= (float
                   (#'sut/γ h obs state-num (nth (sut/states h) state-num)))
                  (float
                   (apply
                    + (map
                       (fn [j] (#'sut/ξ h
                                        obs
                                        state-num
                                        (nth (sut/states h) state-num)
                                        j))
                       (sut/states h))))))
             (range (count (sut/states h)))))))))

(deftest test-baum-welch-1
  (testing "Testing if Baum-Welch converges towards correct hi-prec. value."
    (let [hmm (sut/init-hmm {:s1 1/5 :s2 4/5}
                            {:s1 {:s1 1/2 :s2 1/2}
                             :s2 {:s1 3/10 :s2 7/10}}
                            {:s1 {:N 3/10 :E 7/10}
                             :s2 {:N 4/5 :E 1/5}})]
      (is
       (= 3364526423555802/15688595904862135
          (get-in
           (sut/baum-welch hmm [:N :N :N :N :N :E :E :N :N :N] 1)
           [:A :s2 :s1]))))))

(deftest test-baum-welch-2
  (testing "Testing if Baum-Welch converges towards correct value."
    (let [hmm (sut/init-hmm {:s1 0.2 :s2 0.8}
                            {:s1 {:s1 0.5 :s2 0.5}
                             :s2 {:s1 0.3 :s2 0.7}}
                            {:s1 {:N 0.3 :E 0.7}
                             :s2 {:N 0.8 :E 0.2}})]
      (is
       (= (float 0.14285715)
          (float
           (get-in
            (sut/baum-welch hmm [:N :N :N :N :N :E :E :N :N :N] 1000)
            [:A :s2 :s1])))))))

