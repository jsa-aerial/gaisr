;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                      B I O . A P P S . Q I N                             ;;
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

(ns edu.bc.bio.apps.qin

  "First 'application' level module.  This may not make sense to be
   part of the basic gaisr system.  Should be factored out and then
   use REST calls to gaisr to get the information.

   Investigate the occurance and properties of nc-RNAs in the Qin 2010
   microbiome data sets.  This will search for RNAs (typically
   ribosomal to start), constrained by associated context AND
   discovery of NEW context for these representative RNAs.  New
   context discovery will be driven by SCCS and its Clustering.
  "

  (:require [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.io :as io]
            [incanter.core]
            [incanter.charts]
            [edu.bc.fs :as fs]
            [edu.bc.utils.clustering :as clu])
  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
	edu.bc.bio.sequtils.dists
        ))



(def qin-2ksq-totcnts
     (->> (fs/re-directory-files "/data2/Bio/MetaG1/FastaFiles" "*.seq.fa")
          sort
          (xfold (fn[fa]
                   (->> fa (#(read-seqs % :info :both))
                        (#(freqs&probs
                           1 % (fn[_ es-coll]
                                 (reduce (fn[fm [e s]]
                                           (let [cnt (count s)]
                                             (assoc fm cnt
                                                    (inc (fm cnt 0)))))
                                         {} es-coll))))
                        first (sort-by first >)
                        (take-until (fn[[sz cnt]] (< sz 2000)))
                        (map (fn[[sz cnt]] cnt)) sum)))))

(def qin-totcnts
     (->> (fs/re-directory-files "/data2/Bio/MetaG1/FastaFiles" "*.seq.fa")
          sort
          (xfold (fn[fa]
                   (->> fa (#(read-seqs % :info :both))
                        (#(freqs&probs
                           1 % (fn[_ es-coll]
                                 (reduce (fn[fm [e s]]
                                           (let [cnt (count s)]
                                             (assoc fm cnt
                                                    (inc (fm cnt 0)))))
                                         {} es-coll))))
                        first (sort-by first >)
                        (map (fn[[sz cnt]] cnt)) sum)))))


(def qin-2ksq-totbases
     (->> (fs/re-directory-files "/data2/Bio/MetaG1/FastaFiles" "*.seq.fa")
          sort
          (xfold (fn[fa]
                   (->> fa (#(read-seqs % :info :both))
                        (#(freqs&probs
                           1 % (fn[_ es-coll]
                                 (reduce (fn[fm [e s]]
                                           (let [cnt (count s)]
                                             (assoc fm cnt (inc (fm cnt 0)))))
                                         {} es-coll))))
                        first (sort-by first >)
                        (take-until (fn[[sz cnt]] (< sz 2000)))
                        (sum (fn[[sz cnt]] (* sz cnt))))))))

(def qin-totbases
     (->> (fs/re-directory-files "/data2/Bio/MetaG1/FastaFiles" "*.seq.fa")
          sort
          (xfold (fn[fa]
                   (->> fa (#(read-seqs % :info :both))
                        (#(freqs&probs
                           1 % (fn[_ es-coll]
                                 (reduce (fn[fm [e s]]
                                           (let [cnt (count s)]
                                             (assoc fm cnt (inc (fm cnt 0)))))
                                         {} es-coll))))
                        first (sort-by first >)
                        (sum (fn[[sz cnt]] (* sz cnt))))))))