(ns markov-chains.core-test
  (:require [clj-async-profiler.core :as prof]
            [clojure.test :refer :all]
            [criterium.core :as criterium]
            [markov-chains.core :as sut]
            [no.disassemble :as nda]))

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
;;      (= [{:rain 0.81818181818181810, :no-rain 0.18181818181818182}
;;          {:rain 0.88335704125177800, :no-rain 0.11664295874822190}
;;          {:rain 0.19066793972352525, :no-rain 0.80933206027647480}
;;          {:rain 0.73079400458498200, :no-rain 0.26920599541501794}
;;          {:rain 0.86733888957548470, :no-rain 0.13266111042451528}]
;;         (#'sut/forwards example-hmm-rn
;;                         [:umbrella :umbrella :no-umbrella :umbrella :umbrella])))))

(defn test-fws []
  (criterium.core/report-result
   (criterium.core/quick-benchmark
    (#'sut/forwards (sut/init-hmm {:rain 0.5 :no-rain 0.5}
                                  {:rain {:rain 0.7 :no-rain 0.3}
                                   :no-rain {:rain 0.3 :no-rain 0.7}}
                                  {:rain {:umbrella 0.9 :no-umbrella 0.1}
                                   :no-rain {:umbrella 0.2 :no-umbrella 0.8}})
                    [:umbrella :umbrella :no-umbrella :no-umbrella :umbrella])
    {})))

(defn test-fw []
  (criterium.core/report-result
   (criterium.core/quick-benchmark
    (#'sut/forward (sut/init-hmm {:rain 0.5 :no-rain 0.5}
                                  {:rain {:rain 0.7 :no-rain 0.3}
                                   :no-rain {:rain 0.3 :no-rain 0.7}}
                                  {:rain {:umbrella 0.9 :no-umbrella 0.1}
                                   :no-rain {:umbrella 0.2 :no-umbrella 0.8}})
                   [:umbrella :umbrella :no-umbrella :no-umbrella :umbrella]
                   2
                   :rain)
    {})))

;; Evaluation count : 49446 in 6 samples of 8241 calls.
;;              Execution time mean : 12.524465 µs
;;     Execution time std-deviation : 397.544214 ns
;;    Execution time lower quantile : 12.082332 µs ( 2.5%)
;;    Execution time upper quantile : 13.056886 µs (97.5%)
;;                    Overhead used : 6.376173 ns



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


(defn serve-files []
  (prof/serve-files 8080))
                           
