;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                       U T I L S . T R E E S                              ;;
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

(ns edu.bc.utils.trees

  "Tree branch/node oriented wrapper for zippers."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io])
  (:use edu.bc.utils
        clojure.contrib.math
        [slingshot.slingshot
         :only [try+]]
        [clojure.pprint
         :only [cl-format]]
        ))

;;; *** NOTE: Getting Close, but still NOT TOTALLY COOKED YET!!! ***


(defn map-branch?
  "Predicate defining what constitutes a map branch for map zippers"
  [x]
  (and x (map? x) (not= :!no! (get x :val :!no!))))

(defn map-zip
  "Define a map zipper.  Returns a zipper for nested maps, given a
   starting/root map"
  [root]
  (zip/zipper map-branch?
              #(seq (:children %))
              (fn [node children]
                {:val (:val node) :children (and children (vec children))})
              root))


(defn loc?
  "Return whether x is a valid tree location (aka zipper branch)"
  [x]
  (let [m (meta x)]
    (and m (m :zip/branch?))))

(defn loc-type
  "Return the type of the tree location (zipper branch) loc.  This
   will be one of the extent zipper types (from zip/zipper) available.
   In practice, it is going to be a vector, seq, or map.
  "
  [loc]
  (condp = (loc? loc)
      vector? :vec
      seq?    :list
      map-branch? :map
      nil))

(def ^{:dynamic true
       :doc "Default kind of tree rep - :vec, :list, or :map"}
     *default-tree-kind* :vec)

(defn tree-dispatch-fn [x]
  (if-let [k (loc-type x)]
    k
    (if (keyword? x)
      x
      (cond
       (vector? x) :vec
       (not (coll? x)) *default-tree-kind*
       (map? x) :map
       (seq? x) :list))))


(defn nest-item [loc x valuefn]
  ;;(prn loc x)
  (let [item1 (if (loc? x) (valuefn x) x)]
    (-> loc (zip/append-child item1))))

(defn nest-items [loc nodes valuefn]
  (reduce (fn[loc child] (nest-item loc child valuefn))
          loc
          nodes))

(defn new-loc
  "Return an empty location of type kind an element of
  #{:vec :map :list}"
  ([kind]
     (case kind
           :vec (zip/vector-zip [])
           :map (map-zip {:val nil :children []})
           :list (zip/seq-zip (list))))
  ([kind val]
     (new-loc kind val zip/node))
  ([kind val valuefn]
     (if (= :map kind)
       (zip/edit (new-loc kind)
                 (fn[v] {:val (if (loc? val) (valuefn val) val)
                         :children (v :children)}))
       (nest-item (new-loc kind) val valuefn))))


(defn add-children [loc & children]
  (nest-items loc children zip/node))


(defn tree-branch
  ([x y z & more]
     (let [nodes (->> (seq more) (cons z) (cons y))
           nodes (if (keyword? x) nodes (cons x nodes))
           newloc (new-loc (tree-dispatch-fn x) (first nodes) zip/root)]
       (nest-items newloc (rest nodes) zip/root)))
  ([x y]
     (let [nodes (list y)
           nodes (if (keyword? x) nodes (cons x nodes))
           newloc (new-loc (tree-dispatch-fn x) (first nodes) zip/root)]
       (nest-items newloc (rest nodes) zip/root)))
  ([x]
     (if (keyword? x)
       (new-loc (tree-dispatch-fn x))
       (new-loc (tree-dispatch-fn x) x zip/root))))

(defn tree-node
  ([x y z & more]
     (let [nodes (->> (seq more) (cons z) (cons y))
           nodes (if (keyword? x) nodes (cons x nodes))
           newloc (new-loc (tree-dispatch-fn x) (first nodes) zip/node)]
       (nest-items newloc (rest nodes) zip/node)))
  ([x y]
     (let [nodes (list y)
           nodes (if (keyword? x) nodes (cons x nodes))
           newloc (new-loc (tree-dispatch-fn x) (first nodes) zip/node)]
       (nest-items newloc (rest nodes) zip/node)))
  ([x]
     (if (keyword? x)
       (new-loc (tree-dispatch-fn x))
       (new-loc (tree-dispatch-fn x) x zip/node))))




(defn filter-children [node f]
  (loop [c (zip/down node)
         result []]
    (if (nil? c)
      result
      (recur (zip/right c)
             (if-let [r (f c)] (conj result r) result)))))


(defn tree-find [loc matchfn & {all :all :or {all false}}]
  (loop [loc loc
         result []]
    (cond
     (zip/end? loc)
     (when (seq result) result)

     (matchfn (zip/node loc))
     (if (not all)
       loc
       (recur (zip/next loc)
              (conj result loc)))

     :else
     (recur (zip/next loc) result))))


(defn tree-edit [loc matchfn editfn
                & {all :all zip :zip :or {all false zip true}}]
  (loop [loc loc]
    (cond
     (zip/end? loc)
     (if zip (zip/root loc) loc)

     all
     (recur (zip/next
             (if (matchfn (zip/node loc))
               (zip/edit loc editfn)
               loc)))

     :else
     (if (matchfn (zip/node loc))
       ((if zip zip/root identity) (zip/edit loc editfn))
       (recur (zip/next loc))))))


(defn tree-remove [loc matchfn
                   & {all :all zip :zip
                      :or {all false zip true}}]
  (loop [loc loc]
    (cond
     (zip/end? loc)
     (if zip (zip/root loc) loc)

     all
     (recur (zip/next
             (if (matchfn (zip/node loc))
               (zip/remove loc)
               loc)))

     :else
     (if (matchfn (zip/node loc))
       ((if zip zip/root identity) (zip/remove loc))
       (recur (zip/next loc))))))


(defn tree-transact [tree zipfn & opers]
  (try+
   (loop [ops opers
          loc (zipfn tree)]
     (if (empty? ops)
       (zip/root loc)
       (recur (rest ops)
              ((first ops) loc))))
   (catch #(or (map? %) (set? %)) c
     (prn (:message &throw-context))
     tree)
   (catch Exception e
     (print e)
     tree)))




;;;(tree-edit dz #(= % '*)  :all false)
;;;(-> (tree-node 1 '+  77) (tree-branch '/ dz) zip/root)

(comment
  (def data '[[a * b] + [c * d]])
  (def dz (zip/vector-zip data))

  (def data2 '((a * b) + (c * d)))
  (def d2z (zip/seq-zip data2))
  )

