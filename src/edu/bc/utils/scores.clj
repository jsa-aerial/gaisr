;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                        U T I L S . S C O R E S                           ;;
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

(ns edu.bc.utils.scores

  "Various data test and analysis performance scoring."

  (:require [clojure.contrib.math :as math]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        ))


(defn true-positive-rate
  "Ratio of true positives to the total actual positives, the latter
   being the true positives plus the false negatives (the remaining
   actuals not counted as true): (/ TP (+ TP FN)).  AKA 'hit rate',
   'recall', and 'sensitivity'.
  "
  [tp fn]
  (if (== 0 tp fn)
    (float 1.0) ; no real positives and none predicted => perfect recall
    (float (/ tp (+ tp fn)))))

(defn true-negative-rate
  "Ratio of true negatives to the total actual negatives, the latter
   being the true negatives plus the fase positives (/ TN (+ TN FP)).
   AKA 'specificity'.
   "
  [tn fp]
  (if (== 0 tn fp)
    (float 1.0) ; no real negatives and none predicted
    (float (/ tn (+ tn fp)))))


(defn false-positive-rate
  "Ratio of false positives to the total actual negatives, the latter
   being the true negatives plus the false positives (the remaining
   actuals not counted as false): (/ FP (+ TN FP)) = 1 - 'specificity'

   Where 'specificity' = ratio of true negatives to total actual
   negatives: (/ TN (+ TN FP))
  "
  [fp tn]
  (if (== 0 fp tn)
    (float 0.0)
    (float (/ fp (+ fp tn)))))

(defn false-negative-rate
  "Ratio of false negatives to the total actual positives, the latter
   being the false negatives plus the true positives (the remaining
   actual positives not counted as false): (/ FN (+ FN TP)).  AKA
   'miss rate'.

   *** Somehow this originally named the false POSITIVE rate function
       - how did that happen??
  "
  [fn tp]
  (if (== 0 fn tp)
    (float 0.0)
    (float (/ fn (+ fn tp)))))


(defn acc
  "ACCuracy of classification and (by extension across multiple cases)
   classifier.  TP and TN are the classification TruePostives and
   TrueNegatives respectivel.  AP and AN are the ActualPositives and
   ActualNegatives respectively.  Ratio of sum of correct
   classifications to total values: (TP + TN) / (AP + AN), so result
   lies in [0, 1]
  "
  [tp tn ap an]
  (float (/ (+ tp tn) (+ ap an))))

(defn ppv
  "Positive Predictive Value, aka 'precision'.  Ratio of true
   positives to all predicted positives: TP / (TP + FP).
  "
  [tp fp]
  (if (== 0 tp fp)
    (float 1.0) ; no real positives and none predicted => perfect precision
    (float (/ tp (+ tp fp)))))

(defn mcc
  "Mathews Correlation Coefficient aka 'phi-coefficient'.  Correlation
   indicator (coefficient) between real and predicted classifications.
   Values lie in [-1, 1], where -1 indicates complete negative
   correlation (total disagreement between real and predicted), 0
   indicates prediction no better than random guess, and 1 indicates
   perfect correlation (complete agreement).

   Returns (TP*TN - FP*FN) / sqrt(P*N*P'*N'), where

   P = total real positives = TP + FN
   N = total real negatives = TN + FP
   P' = predicted positives = TP + FP
   N' = predicted negatives = TN + FN
  "
  [tp tn fp fn]
  (let [p  (+ tp fn)
        n  (+ tn fp)
        p' (+ tp fp)
        n' (+ tn fn)
        D  (* p n p' n')]
    (float (/ (- (* tp tn) (* fp fn))
              (if (zero? D) 1 (math/sqrt D))))))



(defn ROC-opt-classifier-pt
  ""
  [pts tprfn fprfn]
  (reduce (fn[[r x] ptrec]
            (let [tpr (tprfn ptrec)
                  fpr (fprfn ptrec)
                  sc (- 1 (- tpr fpr))]
              (if (< sc r)
                [sc ptrec]
                [r x])))
          [10.0 {}] pts))


(defn AUC
  "A simple trapezoidal area under the curve implementation mostly
   intended for ROC curve results.
  "
  [pts]
  (sum (fn[[[xi yi :as A] [xj yj :as B] :as C]]
         (let [[[xi yi] [xj yj]] (if (< yi yj) C [B A])
               base (math/abs (- xj xi))
               height (- yj yi)]
           (+ (* base yi)  ; rectangle
              (* 1/2 base height)))) ; triangle
       pts))
