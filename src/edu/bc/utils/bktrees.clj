;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                       U T I L S . B K T R E E S                          ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011-2012 Trustees of Boston College                       ;;
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

(ns edu.bc.utils.bktrees

  "BK-trees: Metric trees over discrete spaces.  Parameterizable with
   various (not necessarily discrete) metrics and arbitrary values."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io])

  (:use edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.trees

        clojure.contrib.math
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        ))


(defn node [item distance & {val :val :or {val false}}]
  (let [n (tree-node :map {:item item :dist distance})]
    (if val (zip/node n) n)))

(defn item [node]
  (:item (:val (if (loc? node) (zip/node node) node))))

(defn dist [node]
  (:dist (:val (if (loc? node) (zip/node node) node))))

(defn children [node]
  (if (loc? node) (zip/children node) (node :children)))


(defn insert
  ([tree dfn x]
     (->
      (let [d  (dfn x (item tree))]
        (if (= d 0)
          tree
          (if-let [c (first (filter-children tree #(when (= (dist %) d) %)))]
            (insert c dfn x)
            (add-children tree (node x d)))))
      zip/root
      map-zip))

  ([tree dfn x & others]
     (reduce (fn[loc w]
               (insert loc dfn w))
             tree (conj others x))))


(defn query [tree x threshold dfn
             & {:keys [limit ordered] :or {limit 10 ordered true}}]
  (letfn [(bkfind [subtree]
            (let [d (dfn x (item subtree))
                  cr (reduce (fn[v node] ; Reduce children.
                               (if (>= threshold (abs (- d (dist node))))
                                 (concat v (bkfind node))
                                 v))
                             [] (children subtree))]
              (if (<= d threshold)
                (conj cr {:dist d :item (item subtree)})
                cr)))]
    (take limit (sort #(< (:dist %1) (:dist %2))
                      (bkfind (if (loc? tree) tree (map-zip tree)))))))




(defn walk
  "Walk TREE applying F to each node and current result.  RESULT is
   starting value for F and defaults to nil (empty list).  F is a
   function of two parameters: result, node, in that order.  F should
   return a value compatible with result"
  [f tree & [result]]
  (letfn [(_walk [v subtree]
            (let [v (f v subtree)]
              (reduce (fn[v node] (_walk v node))
                      v (children subtree))))]
    (_walk result tree)))


(defn node-count
  "Return the count of nodes in tree"
  [tree]
  (walk (fn[result _] (inc result)) tree 0))


(defn maximum-depth
  "Returns maximum depth of TREE."
  [tree]
  (letfn [(scan [depth subtree]
            (apply max depth
                   (map #(scan (inc depth) %) (children subtree))))]
    (scan 0 tree)))


(defn child-cnts
  "Returns an descending sorted seq of children counts per node of the
   supplied `TREE' as pairs (k cnt), where k is the key of the node
   with cnt children."  [tree]
  (sort #(> (second %1) (second %2))
        (walk (fn[v n]
                (if (seq (children n))
                  (conj v [(item n) (count (children n))])
                  v))
              tree [])))

(defn avg-child-cnt
  "Return average child count per node in TREE."
  [tree]
  (let [ccnts (child-cnts tree)
        numcs (sum (map second ccnts)) ; child count
        numps (count ccnts)]           ; parent count
    (if (zero? numcs) 0 (double (/ numcs numps)))))



;;; -----------------------------------------------------------------------;;;
;;;
(comment
  (def *bktree* (atom {}))

  (let [terms (io/read-lines
               "/other/home/jsa/FG/Mjmy/QuickSite/title-words.txt")]
    (swap! *bktree* (fn[_] (zip/root (node (first terms) 0))))
    (swap! *bktree*
           (fn[v & args]
             (tree-transact
              @*bktree* map-zip
              #(apply insert % levenshtein (rest terms)))))
    :loaded)

  (query @*bktree* "hrassistant" 5 levenshtein :limit 5)
  (query @*bktree* "Contacts" 3 levenshtein :limit 5)

  (->
   (tree-transact
    (zip/root (node "abc" 0))
    map-zip
    #(insert % levenshtein "def" )
    #(insert % levenshtein "happy")
    #(insert % levenshtein "def"))
   node-count)
   (tree-transact
    map-zip
    #(insert % levenshtein "one" "two" "three" "four"))
   node-count)

  )
