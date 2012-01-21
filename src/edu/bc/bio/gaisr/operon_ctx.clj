;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                            O P E R O N - C T X                           ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011 Trustees of Boston College                            ;;
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
;; Authors: Shermin Pei, Jon Anthony                                        ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns edu.bc.bio.gaisr.operon-ctx

  "Operon context matching for sequence 'hits'.  Data sets pulled from
   prebuilt operon tables in extended biosql database"

  (:require [clojure.contrib.sql :as sql]
            [org.bituf.clj-dbcp :as dbcp]
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.contrib.math :as math]
            [clojure.contrib.trace :as tr]
            )

  (:use edu.bc.utils
        [edu.bc.bio.gaisr.db-actions
         :only [operon-location-query]]

        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)])

  (:import javax.sql.DataSource
           com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource))



;;; Two more semantically clear short hands for determining non
;;; symetrical "betweeness" of a number B to two others:
;;;
;;; a < b <= c and a <= b < c
;;;
;;; Maybe these should go into utils??
;;;
(defmacro <_<=
  "a < b <= c"
  [a b c]
  `(let [a# ~a
         b# ~b
         c# ~c]
     (and (< a# b#) (<= b# c#))))

(defmacro <=_< [a b c]
  "a <= b < c"
  `(let [a# ~a
         b# ~b
         c# ~c]
     (and (<= a# b#) (< b# c#))))


;;; minimum and maximum is the min/max over all first or all last
;;; numbers (as per STRAND) in the pairs in the pre sorted vector V
;;;
(defn minimum [strand v]
  (apply min (flatten (if (pos? strand) (first v) (last v)))))

(defn maximum [strand v]
  (apply max (flatten (if (pos? strand) (last v) (first v)))))


;;; x between y and z or y between z and x
;;;
(defn between?
  "x < y <= z or z < y <= x"
  [x y z]
  (or (<_<= x y z) (<=_< z y x)))

;;;picks closest value (y or z) to x
;;;
(defn closer?
  "Returns n that minimizes (abs (- x n)), n in #{y z}"
  [x y z]
  (let [a (min y z)
        b (max y z)]
  (if (<= (math/abs (- x a)) (math/abs(- x b))) a b)))




;;; determines the range to take given an operon
;;;
(defn get-pos [operon from upstream downstream & [p]]
  (loop [i operon
         prev-operon p]
    (let [[cur-start cur-end] (first i)
          next-start (or (first (second i)) cur-end)
          p-end (second prev-operon)]

      ;; (prn :GET-POS :=> from i prev-operon)
      ;; from must be between the end of the end of the previous
      ;; operon and the current operons start (ie inter-operon region)
      ;; unless there is no previous operon
      (cond
       (empty? i)
       (raise :type :bad-operon-seq :args [i prev-operon])

       ;; return [(closer of p-end or from-upstream) (from+25)]
       (between? (or p-end (- from upstream)) from cur-start)
       [(closer? from (- from upstream) p-end) (+ from downstream)]

       ;; else from must be between genes in an operon
       ;;(and (< cur-end from) (<= from next-start))
       (between? cur-end from next-start)
       [cur-end (+ from downstream)]

       ;;from must be in the middle of a gene and the gene is either:
       ;;the first gene in the organism or not
       ;;(and (< cur-start from) (<= from cur-end))
       (between? cur-start from cur-end)
       (if (nil? p-end)
         ;; return range where start must be at least 1
         [(if (pos? (- from upstream)) (- from upstream) 1) (+ from downstream)]
         ;; returns [(closer of p-end or from-upstream) (from+25)]
         [(closer? from (- from upstream) p-end) (+ from downstream)])

       :else
       (recur (rest i)
              (first i))))))


(defn hit-between? [strand min hit-from max]
  (or (and (pos? strand) (<_<= min hit-from max))
      ;;(and (neg? strand) (<_<= max hit-from min))))
      (and (neg? strand) (<=_< min hit-from max))))

(defn find-loc
  ([loc from strand] (find-loc loc from strand 250 25))

  ([loc from strand upstream downstream]
     (when (not (empty? loc))
       (let [sorted-loc (vals (sort-by key loc))
             min-all (minimum strand sorted-loc)
             max-all (maximum strand sorted-loc)]
         ;;(prn sorted-loc)
         ;; hit location must be between the first and last operon of
         ;; the organism or else no operon information can be used
         (cond
          ;; check to make sure that the hit is between the first
          ;; and last operon
          ;; loop through all operons where i = list of all operons remaining
          (hit-between? strand min-all from max-all)
          (loop [i sorted-loc
                 prev-operon []]
            (when (not (empty? i))
              (let [operon (first i)
                    min-opr (minimum strand
                                     (list (concat [prev-operon] operon)))
                    max-opr (maximum strand
                                     (list (concat [prev-operon] operon)))]
                ;; from must be between the end of the last operon and
                ;; the end of current operon
                (if (hit-between? strand min-opr from max-opr)
                  (get-pos operon from upstream downstream prev-operon)
                  (recur (rest i)
                         (last operon))))))

          ;; checks + strand if from is outside the operon range but still
          ;; within (+ max upstream)
          (or (and (pos? strand)
                   (<_<= (maximum strand sorted-loc)
                         from
                         (+ upstream (maximum strand sorted-loc))))
              (and (neg? strand)
                   (<=_< (+ upstream (minimum strand sorted-loc))
                         from
                         (minimum strand sorted-loc))))
          (get-pos [[(- from 10) (+ from 10)]]
                   from upstream downstream
                   (last (last sorted-loc)))

          :else
          ;; no operon information used so a static amount is taken
          ;;(do (prn "no operon info used")
          [(if (pos? (- from upstream)) (- from upstream) 1)
           (if (pos? (+ from downstream)) (+ from downstream) 1)])))))




(defn input [sql-out]
  (map (fn [m] [(m :operon_id) [(m :start) (m :end)] (m :strand)]) sql-out))

(defn add-value [m id loc]
  (conj (get m id []) loc))

(defn make-hash [inseq]
  (reduce (fn [m [id [from to] strand]]
            (let [loc (if (pos? strand) [from to] [to from])]
              (assoc m strand
                     (assoc (get m strand {})
                       id (add-value (get m strand) id loc)))))
          {} inseq))


(defn reverse-map [m]
  (reduce (fn [rmap [id loc]]
            (assoc rmap (* -1 id) (vec (rseq loc))))
            {} m))

(defn get-species [s]
  (operon-location-query s))

(defn species-hash [species]
  (let [h (make-hash (input (get-species species)))]
    (assoc h -1 (reverse-map (get h -1)))))


(defn get-region
  ([species strand hit-start]
     (get-region species strand hit-start 250 25))

  ([species strand hit-start take-upstream take-downstream]
     (let [upstream (* strand take-upstream)
           downstream (* strand take-downstream)]
       (or (find-loc (get (species-hash species) strand)
                     hit-start strand upstream downstream)
           ;;return range where start must be at least 1
           [(if (pos? (- hit-start upstream)) (- hit-start upstream) 1)
            (+ hit-start downstream)]))))



;;; Simple sanity check test.
(defn sanity-check
  ([] (sanity-check "NC_010695"))
  ([species] (sanity-check species 25))
  ([species downstream]
     (let [locs [[1 700 500] [1 1750 250] [1 960 25]
                 [1 8750 250] [1 1500 250]
                 ;; Minus strands
                 [-1 6950 100] [-1 6200 500] [-1 8000 250]
                 [-1 6800 250] [-1 8600 250]]]
       (map (fn [[a b c]]
              (let [result (get-region species a b c downstream)]
                (prn "test result" species [a b c] result)
                [[a b c] :=> result]))
            locs))))

;;;(def abc (get-species "NC_010695"))
;;;(get-region "NC_010695" +1 8229 8000)
;;;
;;; (defn x [strand]
;;;   (get (species-hash "NC_010695") strand))
