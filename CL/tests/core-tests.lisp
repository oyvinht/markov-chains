(asdf:oos 'asdf:load-op 'fiveam)

(in-package :markov-chains)

(defmacro vec (&rest items)
  `(apply #'vector (quote ,items)))

;; Baseline tests for creating a HMM and querying probabilities

#+5am
(defvar *example-state-sequences*
  (list (vec a a a a b c a b d a d b f)
        (vec a a f a b b c d a a c)
        (vec a a b)
        (vec c a c f)
        (vec a a f b d e a)
        (vec f c c a b e d a a f a c)))

#+5am
(defvar *example-state-outcomes*
  (mapcar
   (lambda (el) (apply #'cons el))
   '((a 0)(a 0)(a 0)(c 0)(c 2)(b 0)(b 0)(e 0)(b 0)(a 1)(b 2)(a 1)(c 4)(d 1)(a 1)
     (b 0)(d 6)(e 0)(a 15)(c 4)(c 2)(a 0)(b 8)(a 2)(e 5)(a 2)(d 4)(c 20)(e 9))))

#+5am
(5am:test test-transition-probs
  "Checking probability of transitioning to state b when in a (25%)"
  (let ((hmm (make-hmm *example-state-sequences* *example-state-outcomes*)))
    (5am:is (= (transition-prob hmm 'a 'b) 0.25))))

#+5am
(5am:test test-outcome-probs
  "Checking probability for getting outcome 2 from state b (1/6)"
  (let ((hmm (make-hmm *example-state-sequences* *example-state-outcomes*)))
    (5am:is (= (outcome-prob hmm 'b 2) 1/6))))
