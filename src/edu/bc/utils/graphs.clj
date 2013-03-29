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


(defn dist-matrix
  ""
  [distfn coll & {:keys [sym keyfn] :or {sym true keyfn identity}}]
  (let [cnt (count coll)]
    (reduce
     (fn[M [i j d]]
       (let [ik (keyfn i)
             jk (keyfn j)]
         (if sym (assoc M [ik jk] d [jk ik] d) (assoc M [ik jk] d))))
     {} (xfold (fn [[i j]] [i j (distfn i j)])
               (for [k (range cnt)
                     i (drop k coll)
                     j (if sym (take (- cnt k) coll) coll)
                     :when (not= i j)]
                 [i j])))))


(defn knn [k distfn p coll]
  (take k (sort-by #(distfn p %) coll)))

(defn knn-graph [k distfn coll]
  (let [coll (set coll)]
    (reduce
     (fn[G p] (assoc G p (knn k distfn p (set/difference coll #{p}))))
     {} coll)))


(defn krnn-graph [k distfn coll]
  (let [knngrph (knn-graph k distfn coll)]
    (conj (reduce (fn[[krnnM rnncntM] p]
                    (reduce (fn[[krnnM rnncntM] q]
                              [(assoc krnnM q (conj (get krnnM q []) p))
                               (assoc rnncntM q (inc (get rnncntM q 0)))])
                            [krnnM rnncntM] (knngrph p)))
                  [{} {}] coll)
          knngrph)))


(defn split-krnn [k krnngrph rnncntM knngrph]
  (let [[rnnG>k rnnG<k]
        (reduce (fn[[rnnG>k rnnG<k] p]
                  (if (< (rnncntM p) k)
                    (let [rnnG>k (reduce
                                  (fn[G q]
                                    (if (G q)
                                      (assoc G q (remove #(= % p) (G q)))
                                      G)) ; If it was removed, don't put back!
                                  rnnG>k (knngrph p))]
                      [(dissoc rnnG>k p)
                       (assoc rnnG<k p (rnnG>k p))])
                    [rnnG>k rnnG<k]))
                [krnngrph {}] (keys rnncntM))
        ;; Make final pass on just constructed outliers removing edge
        ;; set points in rnnG>k (well mostly - really any not in
        ;; rnnG<k)
        rnnG<k (let [nodes (keys rnnG<k)]
                 (reduce (fn[rnnG<k q]
                           (if (not (in q nodes))
                             (reduce (fn[G n]
                                       (assoc G n (remove #(= % q) (G n))))
                                     rnnG<k nodes)
                             rnnG<k))
                         rnnG<k (-> rnnG<k vals flatten set)))]
    [rnnG>k rnnG<k]))


(defn refoldin-outliers [krnngrph rnnG>k rnnG<k]
  (let [clusters (first (reduce (fn[[M i] scc]
                                  [(assoc M i scc) (inc i)])
                                [{} 0] rnnG>k))
        outliers (->> rnnG<k (map seq) flatten
                      (map #(do [% (set (krnngrph %))]))
                      (into {}))]
    #_(prn clusters)
    (->>
     (reduce (fn[clus [o ons]]
               (->> clus (map (fn[[k v]] [k (count (set/intersection v ons))]))
                    (sort-by second >) first
                    ((fn[[k cnt]]
                       (if (not= 0 cnt)
                         (assoc clus k (conj (clus k) o))
                         (assoc clus (gen-uid) #{o}))))))
             clusters outliers)
     vals set)))


(defn krnn-clust
  ""
  [k distfn coll & {:keys [keyfn] :or {keyfn identity}}]
  (let [dm (dist-matrix distfn coll :keyfn keyfn)
        kcoll (map keyfn coll)
        [krnngrph rnncntM knngrph] (krnn-graph
                                    k #(get dm [%1 %2]) kcoll)]
    (->> (split-krnn k krnngrph rnncntM knngrph)
         (map #(tarjan (keys %) %))
         (apply refoldin-outliers krnngrph))))






(comment

  (def dm (let [coll (->> "/home/jsa/Bio/FreqDicts/NewRFAM/RF00504-seed-NC.sto"
                          (#(read-seqs % :info :name))
                          (#(get-adjusted-seqs % 1100))
                          (map (fn[i es] [i es]) (iterate inc 0)))
                wz 16
                distfn (fn[[_ [nx sx]] [_ [ny sy]]]
                         (jensen-shannon (probs wz sx) (probs wz sy)))
                keyfn first
                dm (gr/dist-matrix distfn coll :keyfn keyfn)]
            {:coll coll :wz wz :dm dm
             :kcoll (map keyfn coll)
             :distfn distfn :keyfn keyfn}))
(def clu-info
     (let [{:keys [coll wz dm kcoll distfn keyfn]} dm
           distfn2 (fn[l r]
                     (apply jensen-shannon
                            (map #(if (map? %) % (probs wz %))
                                 [l r])))
           avgfn (fn
                   ([sqs]
                      (hybrid-dictionary wz sqs))
                   ([x & xs]
                      (hybrid-dictionary wz (cons x xs))))
           ]
       (for [k (range 4 16)
             :let [[krnngrph rnncntM knngrph]
                   (gr/krnn-graph k #(get dm [%1 %2]) kcoll)]]
         (let [clusters (->> (gr/split-krnn k krnngrph rnncntM knngrph)
                             ;;(#(do (prn :G>k/G<k %) %))
                             (map #(gr/tarjan (keys %) %))
                             ;;(#(do (prn :sccs %) %))
                             (apply gr/refoldin-outliers krnngrph)
                             (#(do (prn :final %) %)))
               ent-clus (let [coll (into {} coll)]
                          (map (fn[scc]
                                 (map (fn[x] (-> x coll first)) scc))
                               clusters))
               seq-clus (let [coll (into {} coll)]
                          (map (fn[scc]
                                 (let [x (xfold (fn[x]
                                                  (->> x coll second
                                                       (probs wz)))
                                                scc)]
                                   [(avgfn x) x]))
                               clusters))]
           [k (clu/S-Dbw-index distfn2 seq-clus :avgfn avgfn) ent-clus]))))



  )


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

  )
