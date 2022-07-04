(in-package :markov-chains)

(defun last-state (path)
  (let ((states (car path))) (elt states (1- (length states)))))

(defun argmax (paths)
  "Return the path having the highest probability value"
  (loop for path in paths with best = (first paths)
        when (> (cdr path) (cdr best))
          do (setf best path)
        finally (return best)))

(defun viterbi (hmm observations)
  "The sequence of states having highest likelihood of explaining observations"
  (flet ((best-augmented-path (paths state)
           (argmax
            (mapcar
             (lambda (p)
               (cons (concatenate 'vector (car p) (vector state))
                     (* (cdr p) (transition-prob hmm (last-state p) state))))
             paths))))
    (loop for o across observations and i from 0
          for paths = (if (= i 0)
                          (mapcar (lambda (state)
                                    (cons (vector state)
                                          (* (init-prob hmm state)
                                             (outcome-prob hmm state o))))
                                  (states hmm))
                          (mapcar
                           (lambda (state)
                             (let ((best (best-augmented-path paths state)))
                               (cons (car best)
                                     (* (cdr best)
                                        (outcome-prob hmm state o)))))
                           (states hmm)))
          finally (return (let ((best (argmax paths)))
                            (values (car best) (cdr best)))))))
