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
            [clojure-csv.core :as csv]
            [incanter.core]
            [incanter.charts]
            [edu.bc.fs :as fs]
            [edu.bc.bio.sequtils.tools :as tools]
            [edu.bc.bio.sequtils.sccs :as sccs])
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


(double (/ (sum qin-2ksq-totcnts) (sum qin-totcnts)))
(double (/ (sum qin-2ksq-totbases) (sum qin-totbases)))


(map #(do [%1 %2])
     (->> (map #(float (/ %1 %2)) qin-2ksq-totcnts qin-totcnts)
          (take 10))
     (->> (map #(float (/ %1 %2)) qin-2ksq-totbases qin-totbases)
          (take 10)))




(def qinmicro-dir "/data2/Bio/QinMicro")


(defn qin-seq-name-first-try [l]
  (->> (subs l 1) (str/split #"_") (take 2) (str/join "_")))

(defn qin-seq-name [l]
  (subs l 1))

(defn qin-entry [l]
  (str/replace-re #"\." "-" l))

;;; "/data2/Bio/MetaG1/FastaFiles/MH0001.seq.fa"
(defn split-qinfna
  [Qinfna & {:keys [base sz] :or {base "/data2/BioData/Qin2010Seqs" sz 2000}}]
  (let [Qindir (fs/join base (->> Qinfna fs/basename (str/split #"\.") first))]
    (when (not (fs/exists? Qindir)) (fs/mkdir Qindir))
    (split-join-fasta-file
     Qinfna :base Qindir :pat #"^>"
     :namefn qin-seq-name
     :entryfn qin-entry
     :testfn (fn[nm sq] (>= (count sq) sz)))))

;;;(doseq [qfna (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")]
;;;  (split-qinfna qfna))


(defn combine-qinfna
  [qin-seq-dir outdir]
  (let [seqfnas (fs/directory-files qin-seq-dir ".fna")
        outfile (fs/join outdir (str (fs/basename qin-seq-dir) ".fna"))]
    (io/with-out-writer outfile
      (doseq [fna seqfnas]
        (doseq [l (io/read-lines fna)]
          (println l))))))

(defn build-qin2k
  [qin-dirdir outdir]
  (let [qin-sq-dirs (fs/directory-files qin-dirdir "")]
    (when (not (fs/exists? outdir)) (fs/mkdir outdir))
    (doseq [qdir qin-sq-dirs]
      (combine-qinfna qdir outdir))))

;;;(build-qin2k "/data2/BioData/Qin2010Seqs" (fs/join qinmicro-dir "Qin2K"))


(defn qin-genome-fasta-dir [nm]
  (let [base "/data2/BioData/Qin2010Seqs"
        nm (first (entry-parts (str/replace-re #"\." "-" nm)))
        [x y fname] (str/split #"_" nm)]
    (fs/join base fname)))

(defn ref58-qin-genome-fasta-dir [nm]
  (if (re-find #"^N" nm)
    (refseq58-genome-fasta-dir nm)
    (qin-genome-fasta-dir nm)))

(swap! genome-db-dir-map
       #(assoc % :qin2010 qin-genome-fasta-dir))
(swap! genome-db-dir-map
       #(assoc % :refseq58-qin2010 ref58-qin-genome-fasta-dir))


(defn build-config-files
  []
  (let [dirs (->> (fs/directory-files qinmicro-dir "")
                  (filter #(not (or (re-find #"may" %)
                                    (re-find #"2K" %)
                                    (re-find #"run" %))))
                  sort)
        runtxts (->> (fs/directory-files qinmicro-dir ".txt")
                     (map #(do [(fs/basename %) (io/read-lines %)])))]
    (doseq [d dirs]
      (let [dnm (fs/basename d)
            xpart (str/lower-case dnm)
            sto (->> d (#(fs/directory-files % ".sto"))
                     (remove #(re-find #"(IBD|UC|CD|MH)" %))
                     first fs/basename
                     (str/replace-re #"\.sto" ""))]
        (doseq [[nm lines] runtxts]
          (let [nm (fs/join d (str/replace-re #"xxx" xpart nm))]
            (io/with-out-writer nm
              (doseq [l lines]
                (->> l
                     (str/replace-re #"RNAXXX" sto)
                     (str/replace-re #"XXX" dnm)
                     println)))))))))
(build-config-files)


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


(defn false-negative-rate
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

(defn acc
  "ACCuracy of classification and (by extension across multiple cases)
   classifier.  TP and TN are the classification TruePostives and
   TrueNegatives respectivel.  AP and AN are the ActualPositives and
   ActualNegatives respectively.  Ratio of sum of correct
   classifications to total values: (TP + TN) / (P + N), so result
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
              (if (zero? D) 1 (sqrt D))))))



(def eval-cut-points
     [1.0e-23, 1.0e-21, 1.0e-19, 1.0e-15, 1.0e-11,
      1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 0.001, 0.01
      0.1 0.2 0.3 0.5 0.7 0.9 1.0 2.0 3.0 4.0 5.0 7.0])

(count eval-cut-points)

(defn build-eval-pos-neg-sets
  [csv-hit-file]
  (let [evset (->> csv-hit-file
                   get-csv-entry-info
                   (map (fn[[n s e ev & _]]
                          (let [[s e sd] (if (> (Integer. s) (Integer. e))
                                           [e s -1]
                                           [s e 1])]
                            [(make-entry n s e sd) (Float. ev)])))
                   (sort-by second))
        basedir (->> csv-hit-file
                     fs/split
                     (take-until #(= % "CSV"))
                     ((fn[bits] (apply fs/join (conj (vec bits) "EVAL-ROC")))))]
    (when (not (fs/exists? basedir)) (fs/mkdirs basedir))
    (doseq [cpt eval-cut-points]
      (gen-entry-file
       (map first (take-until (fn[[n ev]] (> ev cpt)) evset))
       (fs/join basedir (str "ev-" cpt "-pos.ent"))))
    (doseq [cpt eval-cut-points]
      (gen-entry-file
       (map first (drop-until (fn[[n ev]] (> ev cpt)) evset))
       (fs/join basedir (str "ev-" cpt "-neg.ent"))))))

(do
  (build-eval-pos-neg-sets
   "/home/kaila/Bio/Test/ROC-data/L20/CSV/EV10/L20-auto-0.sto.all-custom-db.cmsearch.csv")
  (build-eval-pos-neg-sets
   "/home/kaila/Bio/Test/ROC-data/S4/CSV/EV10/S4-auto-0.sto.all-custom-db.cmsearch.csv")
  (build-eval-pos-neg-sets
   "/home/kaila/Bio/Test/ROC-data/S15/CSV/EV10/S15-auto-0.sto.all-custom-db.cmsearch.csv")
  (build-eval-pos-neg-sets
 "/home/kaila/Bio/Test/ROC-data/S6S/CSV/Start/S6-start.sto.all-custom-db.cmsearch.csv")
  (build-eval-pos-neg-sets
 "/home/kaila/Bio/Test/ROC-data/S6M/CSV/Mid/S6-mid.sto.aggregate.cmsearch.csv"))




(defn eval-roc-data
  [& {:keys [base RNAs]
      :or {base "/home/kaila/Bio/Test/ROC-data"
           RNAs ["L20" "S15" "S4" "S6S" "S6M"]}}]
  (let [cutpts eval-cut-points]
    (for [rna RNAs]
      (let [totdir (fs/join base rna "TotPosNeg")
            tot-pos (fs/join totdir (str rna "-pos.ent"))
            tot-neg (fs/join totdir (str rna "-neg.ent"))
            P (-> tot-pos get-entries count)
            N (-> tot-neg get-entries count)]
        (for [cpt cutpts
              :let [dir (fs/join base rna "EVAL-ROC")
                    pos-ent (fs/join dir (str "ev-" cpt "-pos.ent"))
                    neg-ent (fs/join dir (str "ev-" cpt "-neg.ent"))]]
          (let [tp (count (tools/entry-file-intersect true pos-ent tot-pos))
                fn (count (tools/entry-file-intersect true neg-ent tot-pos))
                fp (count (tools/entry-file-intersect true pos-ent tot-neg))
                tn (count (tools/entry-file-intersect true neg-ent tot-neg))]
            [rna cpt
             {:TPR (true-positive-rate tp fn)
              :FPR (false-negative-rate fp tn)
              :MCC (mcc tp tn fp fn)
              :ACC (acc tp tn P N)
              :ppv (ppv tp fp)}])) ))))

(eval-roc-data)




(defparameter sccs-cut-points
     [0.1, 0.2, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
      0.9, 0.925 0.93 0.935 0.94 0.945 0.95 0.955 0.96 0.97
      0.975, 0.98, 0.985, 0.99, 0.995])

(count sccs-cut-points)


(defn perform-sccs
  [csv-hit-file Mre bits basedir outdir csvdir chart-dir]
  (let [bits (fs/split csv-hit-file)
        basedir (->> bits (take-until #(= % "CSV")) (apply fs/join))
        clusters (first (fs/glob (fs/join basedir "CLU-*")))
        stos (ensure-vec
              (if clusters
                (sort (fs/directory-files clusters ".sto"))
                (fs/join basedir
                         (->> bits last
                              (str/split #"\.") first
                              (#(str % ".sto"))))))]
    (doseq [sto stos]
      (sccs/compute-candidate-sets
       sto csv-hit-file
       1 (sccs/get-saved-ctx-size sto)
       :refn jensen-shannon
       :xlate +RY-XLATE+ :alpha ["R" "Y"]
       :refdb :refseq58 :stodb :refseq58
       :crecut 0.01 :limit 19
       :Dy 0.49 :Mre Mre
       :plot-dists chart-dir))

    (when clusters
      (let [stonm (->> bits last (str/split #"\.") first)
            aggr (str stonm "-aggr")
            _ (sccs/aggregate-sccs-ents csvdir csv-hit-file csvdir aggr)
            candidate-clu-ents  (fs/glob (fs/join  csvdir "*final*.ent"))
            hitonly-clu-ents (fs/glob (fs/join  csvdir "*hitonly*.ent"))]
        (doseq [f (concat candidate-clu-ents hitonly-clu-ents)] (fs/rm f))))

    (let [[neg pos] (sort
                     (if clusters
                       (fs/re-directory-files csvdir #"(neg|pos)\.ent$")
                       (fs/glob (fs/join csvdir "*final*.ent"))))
          pos-nm (fs/join outdir (str "sccs-" Mre "-pos.ent"))
          neg-nm (fs/join outdir (str "sccs-" Mre "-neg.ent"))]
      (fs/rename pos pos-nm)
      (fs/rename neg neg-nm))))


(defn build-sccs-pos-neg-sets
  [csv-hit-file]
  (let [bits (fs/split csv-hit-file)
        basedir (->> bits (take-until #(= % "CSV")) (apply fs/join))
        outdir (->> bits
                    (take-until #(= % "CSV"))
                    ((fn[bits] (apply fs/join (conj (vec bits) "SCCS-ROC")))))
        csvdir (apply fs/join basedir
                      (->> bits
                           (drop-until #(= % "CSV"))
                           (take 2)))
        chart-dir (fs/join basedir "Charts")]
    (when (not (fs/exists? outdir)) (fs/mkdirs outdir))
    (when (not (fs/exists? chart-dir)) (fs/mkdirs chart-dir))
    (doseq [Mre sccs-cut-points]
      (println :Mre Mre)
      (perform-sccs csv-hit-file Mre bits basedir outdir csvdir chart-dir)
      (when (> (fd-use) 9000) ; BOGUS??  Why are we leaking FDs???
        (force-gc-finalize)))
    (force-gc-finalize)))

;;;; JVM arg -XX:-MaxFDLimit to allow going to os/proc limit of fds


(do
  (build-sccs-pos-neg-sets
   "/home/kaila/Bio/Test/ROC-data/L20/CSV/EV10/L20-auto-0.sto.all-custom-db.cmsearch.csv" "EV10")
  (build-sccs-pos-neg-sets
   "/home/kaila/Bio/Test/ROC-data/S4/CSV/EV10/S4-auto-0.sto.all-custom-db.cmsearch.csv" "EV10")
  (build-sccs-pos-neg-sets
   "/home/kaila/Bio/Test/ROC-data/S15/CSV/EV10/S15-auto-0.sto.all-custom-db.cmsearch.csv" "EV10")
  (build-sccs-pos-neg-sets "/home/kaila/Bio/Test/ROC-data/S6S/CSV/Start/S6-start.sto.all-custom-db.cmsearch.csv")
  (build-sccs-pos-neg-sets "/home/kaila/Bio/Test/ROC-data/S6M/CSV/Mid/S6-mid.sto.aggregate.cmsearch.csv"))




(defn sccs-roc-data [& {:keys [base RNAs]
                        :or {base "/home/kaila/Bio/Test/ROC-data"
                             RNAs ["L20" "S15" "S4" "S6S" "S6M"]}}]
  (let [cutpts sccs-cut-points]
    (for [rna RNAs]
      (let [totdir (fs/join base rna "TotPosNeg")
            tot-pos (fs/join totdir (str rna "-pos.ent"))
            tot-neg (fs/join totdir (str rna "-neg.ent"))
            P (-> tot-pos get-entries count)
            N (-> tot-neg get-entries count)]
        (for [cpt cutpts
              :let [dir (fs/join base rna "SCCS-ROC")
                    pos-ent (fs/join dir (str "sccs-" cpt "-pos.ent"))
                    neg-ent (fs/join dir (str "sccs-" cpt "-neg.ent"))]]
          (let [tp (count (tools/entry-file-intersect true pos-ent tot-pos))
                fn (count (tools/entry-file-intersect true neg-ent tot-pos))
                fp (count (tools/entry-file-intersect true pos-ent tot-neg))
                tn (count (tools/entry-file-intersect true neg-ent tot-neg))]
            [rna cpt
             {:TPR (true-positive-rate tp fn)
              :FPR (false-negative-rate fp tn)
              :MCC (mcc tp tn fp fn)
              :ACC (acc tp tn P N)
              :ppv (ppv tp fp)}])) ))))


(defn build-chart [rna sccs-data eval-data]
  (let [sccs-data (map last sccs-data)
        eval-data (map last eval-data)
        [sccs-x sccs-y] [(map :FPR sccs-data) (map :TPR sccs-data)]
        [eval-x eval-y] [(map :FPR eval-data) (map :TPR eval-data)]
        [x y] [(range 0.0 1.0 0.1) (range 0.0 1.0 0.1)]
        chart (incanter.charts/scatter-plot
               sccs-x sccs-y
               :x-label "FPR (1-specificity)"
               :y-label "TPR (sensitivity)"
               :title (str rna " SCCS and Eval ROC plots")
               :series-label "SCCS" :legend true)
        chart (incanter.charts/add-points
               chart eval-x eval-y :series-label "Eval")
        chart (incanter.charts/add-points
               chart x y :series-label "Random Guess")]
    chart))

(defn build-roc-curves
  [& {:keys [base RNAs]
      :or {base "/home/kaila/Bio/Test/ROC-data"
           RNAs ["L20" "S15" "S4" "S6S" "S6M"]}}]
  (let [sccs-data (sccs-roc-data :base base :RNAs RNAs)
        eval-data (eval-roc-data :base base :RNAs RNAs)]
    (doseq [i (range 0 (count RNAs))]
      (let [chart (build-chart (nth RNAs i) (nth sccs-data i) (nth eval-data i))
            chart-file (fs/join base "Plots"
                                (str (nth RNAs i) "-roc-plot.png"))]
        (incanter.core/view chart)
        (incanter.core/save
         chart chart-file :width 500 :height 500)))))

(build-roc-curves)


(defn csv-roc-data
  [& {:keys [base RNAs]
      :or {base "/home/kaila/Bio/Test/ROC-data"
           RNAs ["L20" "S15" "S4" "S6S" "S6M"]}}]
  (let [sccs-data (sccs-roc-data :base base :RNAs RNAs)
        eval-data (eval-roc-data :base base :RNAs RNAs)
        step (/ 1.0 (->> sccs-data first count))
        cols ["RNA" "Cutpt" "TPR" "FPR" "Type"]
        guess (map #(do ["" "" (str %) (str %) "Guess"]) (range 0.0 1.0 step))
        xf #(map (fn[[nm cpt vm]]
                   (cons nm (map str [cpt (vm :TPR) (vm :FPR) %1])))
                 %2)]
    (println step guess)
    (doseq [i (range 0 (count RNAs))]
      (let [sccsi (xf "SCCS" (nth sccs-data i))
            evali (xf "EVAL" (nth eval-data i))
            all (concat sccsi evali guess)]
        (spit (fs/join base (str (ffirst all) "-both.csv"))
              (csv/write-csv (cons cols all)))))))

(csv-roc-data :RNAs ["S6S" "S6M"])

(count (map #(do ["" "" % % "Guess"]) (range 0.0 1.0 0.05)))

(map #(do ["" "" (str %) (str %) "Guess"]) (range 0.0 1.0 0.05))





;;; ------------------- Statistical RNA Search ----------------------------;;;


(def bp-scores
     {[\G \C] 1.00, "GC" 1.00, [\C \G] 1.00, "CG" 1.00,
      [\A \U] 0.95, "AU" 0.95, [\U \A] 0.95, "UA" 0.95,
      [\G \U] 0.80, "GU" 0.80, [\U \G] 0.80, "UG" 0.80
      })

(defn base-candidates
  [[cnt base-index-pairs] & {:keys [min-dist] :or {min-dist 4}}]
  (keep (fn[[[b1 i] [b2 j]]]
            (let [bp-score (bp-scores [b1 b2])
                  d (inc (abs (- j i)))]
              (when (and bp-score (> d min-dist))
                (let [n (/ d cnt)
                      score (* n bp-score)]
                [[b1 i] [b2 j] n score]))))
        (combins 2 base-index-pairs)))


(defn plot-score-dist
  [chart-file score-count-pairs]
  (let [[xs ys] [(map first score-count-pairs) (map second score-count-pairs)]
        chart (incanter.charts/scatter-plot
               xs ys
               :x-label "Base pairing score"
               :y-label "count"
               :title "Base Pair Scoring Distribution")]
    (incanter.core/save
     chart chart-file :width 500 :height 500)))



(->> "/data2/Bio/QinMicro/RF00059/RF00059-seed-NC.sto"
     read-seqs
     degap-seqs
     norm-elements
     ;(take 5)
     (map (fn[sq] [(count sq) (map #(do [%1 %2]) sq (iterate inc 0))]))
     (map base-candidates)
     (map (fn[xs] (sort-by #(-> % first second) xs)))
     (map (fn[xs] (map #(/ (round (* % 1.0e15)) 1.0e15) (map last xs))))
     (map (fn[xs] (freqn 1 xs)))
     (map (fn[dist] (sort-by key dist)))
     (map (fn[i dist]
            (plot-score-dist
             (str "/data2/Bio/BPDistRNA/test-bp-dist-" (str i ".png"))
             dist))
          (iterate inc 1)))


;;; 'Control' - random seq behavior
(->> (map #(first (gen-random-bioseq :rna %)) (repeat 10 108))
     ;(take 5)
     (map (fn[sq] [(count sq) (map #(do [%1 %2]) sq (iterate inc 0))]))
     (map base-candidates)
     (map (fn[xs] (sort-by #(-> % first second) xs)))
     (map (fn[xs] (map #(/ (round (* % 1.0e15)) 1.0e15) (map last xs))))
     (map (fn[xs] (freqn 1 xs)))
     (map (fn[dist] (sort-by key dist)))
     (map (fn[i dist]
            (plot-score-dist
             (str "/data2/Bio/BPDistRNA/RAND-bp-dist-" (str i ".png"))
             dist))
          (iterate inc 1)))




(/ (round 123.534545) 1000.0)








(comment

(defn eval-pos-neg-sets
  [csv-hit-file cut-point]
  (let [evset (->> csv-hit-file
                   get-csv-entry-info
                   (map (fn[[n s e ev & _]]
                          (let [[s e sd] (if (> (Integer. s) (Integer. e))
                                           [e s -1]
                                           [s e 1])]
                            [(make-entry n s e sd) (Float. ev)])))
                   (sort-by second))]
    [(take-until (fn[[n ev]] (> ev cut-point)) evset)
     (drop-until (fn[[n ev]] (> ev cut-point)) evset)]))


(defn sccs-pos-neg-sets
  ([sccs-pos-file sccs-neg-file cut-point]
     (let [pos (->> sccs-pos-file
                    slurp csv/parse-csv butlast
                    (map (fn[[nm re]] [nm (Float. re)]))
                    set)
           neg (->> sccs-neg-file
                    slurp csv/parse-csv butlast
                    (map (fn[[nm re]] [nm (Float. re)]))
                    set)
           all (sort-by second (set/union pos neg))]
       [(take-until (fn[[n re]] (> re cut-point)) all)
        (drop-until (fn[[n re]] (> re cut-point)) all)]))
  ([hit-dir cut-point]
     (let [pos (first (fs/glob (fs/join hit-dir "*final.ent")))
           neg (first (fs/glob (fs/join hit-dir "*final-bad.ent")))]
       (sccs-pos-neg-sets pos neg cut-point))))

)