;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                       U T I L S . G R A P H S                            ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011-2013 Trustees of Boston College                       ;;
;;                                                                          ;;
;; Permission is hereby granted, free of charge, to any person obtaining    ;;
;; a copy of this software and associated documentation files (the          ;;
;; "Software"), to deal in the Software without restriction, including      ;;
;; without limitation the rights to use, copy, modify, merge, publish,      ;;
;; distribute, sublicense, and/or sell copies of the Software, and to       ;;
;; permit persons to whom the Software is furnished to do so, subject to    ;;
;; the following conditions:                                                ;;
;;                                                                          ;;
;; The above copyright notice and this permission notice shall be           ;;
;; included in all copies or substantial portions of the Software.          ;;
;;                                                                          ;;
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,          ;;
;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF       ;;
;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                    ;;
;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE   ;;
;; LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ;;
;; OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION    ;;
;; WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          ;;
;;                                                                          ;;
;; Author: Jon Anthony                                                      ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns edu.bc.utils.graphs

  "Various graph algorithms, techniques, functions.
   Generally applicable, but typically used for sequence similarity,
   clustering, and path metrics of various sorts."

  (:require [clojure.contrib.math :as math]
            [clojure.contrib.string :as str]
            [clojure.contrib.graph :as gr]
            [clojure.set :as set]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        edu.bc.utils.probs-stats
        [clojure.pprint
         :only [cl-format]]))



(defn tarjan
  "Returns the strongly connected components of a graph specified by
   its nodes and a successor function succs from node to nodes.  The
   used algorithm is Tarjan's one.
  "
  [nodes succs]
  (letfn [(sc [env node]
              ;; env is a map from nodes to stack length or nil,
              ;; nil means the node is known to belong to another SCC
              ;; there are two special keys: ::stack for the current stack
              ;; and ::sccs for the current set of SCCs
              (if (contains? env node)
                env
                (let [stack (::stack env)
                      n (count stack)
                      env (assoc env node n ::stack (conj stack node))
                      env (reduce (fn [env succ]
                                    (let [env (sc env succ)]
                                      (assoc env node
                                             (min (or (env succ) n)
                                                  (env node)))))
                                  env (succs node))]
                  (if (= n (env node))
                    ;; no link below us in the stack, call it a SCC
                    (let [nodes (::stack env)
                          scc (set (take (- (count nodes) n) nodes))
                          ;; clear all stack lengths for these nodes
                          ;; since this SCC is done
                          env (reduce #(assoc %1 %2 nil) env scc)]
                      (assoc env ::stack stack ::sccs (conj (::sccs env) scc)))
                    env))))]
    (::sccs (reduce sc {::stack () ::sccs #{}} nodes))))




(defn tree->prufer
  "Convert a labelled tree to a Prüfer Sequence."
  ([graph] (tree->prufer graph []))
  ([graph prufer-seq]
     (loop [g graph
            psq prufer-seq]
       (if (= (count g) 2)
         (reverse psq)
         (let [[k v] (apply min-key
                            first (filter (fn[[k v]] (= 1 (count v))) g))]
           (recur (reduce (fn[g [x v]] (assoc g x (disj v k)))
                          {} (dissoc g k))
                  (cons (first v) psq)))))))

(defn prufer->tree
  ""
  [prufer-seq]
  (loop [s (vec prufer-seq)
         l (vec (range 1 (+ (count prufer-seq) 2 1)))
         g (into {} (map #(do [% #{}]) l))]
    (if (= 2 (count l))
      (let [[x y] l]
        (assoc g x (conj (g x #{}) y) y (conj (g y #{}) x)))
      (let [x (first s)
            y (loop [l l]
                (let [n (first l)] (if (not (in n s)) n (recur (rest l)))))]
        (recur (rest s)
               (remove #(= % y) l)
               (assoc g x (conj (g x) y) y (conj (g y) x)))))))






(comment

  (def test-graph
       {:a #{:b}
        :b #{:c}
        :c #{:a}
        :d #{:b :c :e}
        :e #{:d :f}
        :f #{:c :g}
        :g #{:f}
        :h #{:e :g :h}})

  (tarjan (keys test-graph) test-graph)

  (def test-ccgraph
       (struct gr/directed-graph (keys test-graph) test-graph))

  (gr/scc test-ccgraph)

  (def ps-test-g
       {1 #{8 2 7}, 2 #{1}, 3 #{7}, 4 #{5},
        5 #{7 4}, 6 #{7}, 7 #{1 3 6 5}, 8 #{1}})

  (prufer->tree (tree->prufer ps-test-g))
  (reduce (fn[sm n] (assoc sm n (inc (sm n 0))))
          {} (tree->prufer ps-test-g))
  )
