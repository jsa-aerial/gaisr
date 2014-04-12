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
            [edu.bc.bio.sequtils.sccs :as sccs]
            [edu.bc.bio.gaisr.db-actions :as db])
  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.scores
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.dists
        ))

(->> "/home/kaila/Bio/Test/ROC-data/Ecoli/S15/Charts/S15-start-ctxsz.csv"
     slurp csv/parse-csv butlast rest
     (map (fn[[lstg restg nm]] [(Integer. lstg) (Float. restg)]))
     (map (fn[[l re]] [l (expt re (-> l log))])) (sort-by second) (take 10))

(->> "/home/kaila/Bio/Test/ROC-data/Ecoli/S8/Charts/S8-start-ctxsz.csv"
     slurp csv/parse-csv butlast rest
     (map (fn[[lstg restg nm]] [(Integer. lstg) (Float. restg)]))
     (map (fn[[l re]] [l (expt re (-> l log))])) (sort-by second) (take 20))


(->> "/home/kaila/Bio/Test/ROC-data/Sim//RF00174/CLU/Charts/clu-k41-1-ctxsz.csv"
     slurp csv/parse-csv butlast rest
     (map (fn[[lstg restg nm]] [(Integer. lstg) (Float. restg)]))
     (map (fn[[l re]] [l (expt re (-> l log))])) (sort-by second) (take 20))


(->> "/home/kaila/Bio/Test/ROC-data/Sim/RF00379/CLU/Charts/clu-k10-3-ctxsz.csv"
     slurp csv/parse-csv butlast rest
     (map (fn[[lstg restg nm]] [(Integer. lstg) (Float. restg)]))
     (map second)
     ((fn[coll] (let [mn (apply min coll)
                      mx (apply max coll)]
                  [mn mx (- mx mn) (mean coll)]))))

(let [base "/home/kaila/Bio/Test/ROC-data/Sim//RF00634/CLU/Charts"]
  (for [csv (-> base (fs/join "*ctxsz.csv") fs/glob)]
    (->> csv
         slurp csv/parse-csv butlast rest
         (map (fn[[lstg restg nm]] [(Integer. lstg) (Float. restg)]))
         (map (fn[[l re]] [l (expt re (-> l log))]))
         (sort-by second) (take 20))))


(map #(map (fn[ent]
             (let [[nm [s e] st] (entry-parts ent)]
               (make-entry nm (- s 50) (+ e 100) st)))))

(with-genome-db (@genome-db-dir-map :refseq58)
  (->> "/home/kaila/Bio/Test/ROC-data/Sim/RF*" fs/glob sort
       (map #(first (fs/glob (fs/join % "*NC.sto"))))
       (map #(do [% (read-seqs % :info :name)]))
       (map (fn[[f ents]]
              [f (gen-name-seq-pairs ents :ldelta 50 :rdelta 100)]))
       (map (fn[[f pairs]]
              (nms-sqs->fasta-file pairs (fs/replace-type f ".fna"))))
       (#(fs/move % "/home/kaila/Bio/Test/ROC-data/Sim/StartFastas"))))


;;gaisr entry-file-set intersect -fuzzy 103 CSV/RF00050-seed-NC.sto.all-sim-db.cmsearch.csv RF00050-seed-NC.sto

;;gaisr entry-file-set difference ../CSV/RF00050-seed-NC.sto.all-sim-db.cmsearch.csv RF00050-pos.ent

;;; In RF00162, NC_003869/1750321-1750369/-1 is really a duplicate of
;;; NC_003869/1750270-1750369/-1 in the SEED STO!!  So, this will
;;; cause at least one false positive.
;;;
;;; In RF00168, there is another dup
;;;
;;; In RF00504, there are seven missing/dups
;;;
;;; In RF00557 missing or likely dup of 1
;;;
;;;

;;;RF00114/CLU-G:
;;;RF00521/CLU-G:
;;;RF00555/CLU-G:



(defn sto->ent [stofile & {:keys [entfile]}]
  (let [entfile (if entfile
                  (fs/fullpath entfile)
                  (fs/replace-type stofile ".ent"))]
    (io/with-out-writer entfile
      (doseq [ent (read-seqs stofile :info :name)]
        (println ent)))))


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

(keys @genome-db-dir-map)

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


(defparameter eval-cut-points
  [1.0e-33, 1.0e-30, 1.0e-29 1.0e-28 1.0e-27, 1.0e-24, 1.0e-21, 1.0e-18,
   1.0e-15, 1.0e-11, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4,
   0.001, 0.01 0.1 0.2 0.3 0.5 0.7 0.9 1.0])

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
    (when (fs/exists? basedir) (fs/rm-rf basedir))
    (fs/mkdirs basedir)
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

;;; All Ecoli
(doseq [f (fs/glob "/home/kaila/Bio/Test/ROC-data/Ecoli/S*")]
  (let [csv (-> f (fs/join "CSV/*cmsearch.csv") fs/glob first)]
    (build-eval-pos-neg-sets csv)))

;;; All Sim
(doseq [f (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")]
  (let [csv (-> f (fs/join "CSV/*cmsearch.csv") fs/glob first)]
    (build-eval-pos-neg-sets csv)))



(defn eval-roc-data
  [& {:keys [base RNAs]
      :or {base "/home/kaila/Bio/Test/ROC-data"
           RNAs ["L20" "S15" "S4" "S6S" "S6M"]}}]
  (for [rna RNAs]
    (let [dir (fs/join base rna "EVAL-ROC")
          cutpts (->> (fs/glob (fs/join dir "/ev*-pos.ent"))
                      (map fs/basename)
                      (map #(->> % (str/split #"-")
                                 rest butlast
                                 (str/join "-"))))
          totdir (fs/join base rna "TotPosNeg")
          tot-pos (fs/join totdir (str rna "-pos.ent"))
          tot-neg (fs/join totdir (str rna "-neg.ent"))
          P (-> tot-pos get-entries count)
          N (-> tot-neg get-entries count)]
      (for [cpt cutpts
            :let [pos-ent (fs/join dir (str "ev-" cpt "-pos.ent"))
                  neg-ent (fs/join dir (str "ev-" cpt "-neg.ent"))]]
        (let [tp (count (tools/entry-file-intersect true pos-ent tot-pos))
              fn (count (tools/entry-file-intersect true neg-ent tot-pos))
              fp (count (tools/entry-file-intersect true pos-ent tot-neg))
              tn (count (tools/entry-file-intersect true neg-ent tot-neg))]
          [rna cpt
           {:TPR (true-positive-rate tp fn)
            :FPR (false-positive-rate fp tn)
            :MCC (mcc tp tn fp fn)
            :ACC (acc tp tn P N)
            :ppv (ppv tp fp)}])) )))

(eval-roc-data)
(def eval-ecoli
     (eval-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Ecoli"
      :RNAs ["S1" "S4" "S7" "S8" "S15"]))

(->> (range 0.25 0.55 0.01) (map #(do [% (logistic %)])) (take 25))
(->> (range 1.9 3.0 0.04) (map #(do [% (logistic %)])) (take 25))
(->> (range 1.8 3.0 0.03) (map #(do [% (logistic %)])) (take 25))
(->> (range 2.42 3.0 0.03) (map logistic) (take 6))

(defparameter sccs-cut-points
  #_(->> (range -3.7 6 0.4) (map logistic))
  (->> (range -0.2 6 0.2) (map logistic) (take 25)) ; single
  #_(->> (range -0.25 4.0 0.1) (map logistic) (take 25)) ; multi
  #_[0.1, 0.2, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, ; by hand
     0.9, 0.925 0.93 0.935 0.94 0.945 0.95 0.955 0.96 0.97
     0.975, 0.98, 0.985, 0.99, 0.995])

(->> (range -3.7 6 0.4) (map logistic)) count)

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
       :plot-dists chart-dir)
      (when (> (fd-use) 9000) ; BOGUS??  Why are we leaking FDs???
        (force-gc-finalize)))

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
                           (take 1))) ;; *** HACK
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

(build-sccs-pos-neg-sets (->> (fs/directory-files "/home/kaila/Bio/Test/ROC-data/Sim/RF00050/CSV" ".cmsearch.csv") first))

;;;Single Sim
(binding [sccs-cut-points [0.9370266439430035]]
  (time (doall (build-sccs-pos-neg-sets
                (->> (fs/directory-files
                      "/home/kaila/Bio/Test/ROC-data/Sim/RF00559/CSV"
                      ".cmsearch.csv") first)))))

;;; Single Ecoli
(binding [sccs-cut-points (->> (range -0.2 6 0.2) (map logistic) (take 25))]
  (let [base "/home/kaila/Bio/Test/ROC-data/Ecoli"]
    (doseq [rf ["S1"]]
      (doall (build-sccs-pos-neg-sets
              (->> (fs/directory-files
                    (fs/join base rf "CSV" ) ".cmsearch.csv")
                   first))))))

;;; All Ecoli
(binding [sccs-cut-points (->> (range -0.2 6 0.2) (map logistic) (take 25))]
  (let [base "/home/kaila/Bio/Test/ROC-data/Ecoli"]
    (doseq [rf (sort (->> "/home/kaila/Bio/Test/ROC-data/Ecoli/S*"
                          fs/glob (map fs/basename)))]
      (when (fs/directory? (fs/join base rf "SCCS-ROC"))
        (fs/rename (fs/join base rf "SCCS-ROC")
                   (fs/join base rf "ORIG-SCCS-ROC")))
      (time (doall (build-sccs-pos-neg-sets
                    (->> (fs/directory-files
                          (fs/join base rf "CSV" ) ".cmsearch.csv")
                         first)))))))

;;; All Multi Context Sim
(binding [sccs-cut-points (->> (range -0.25 4.0 0.1) (map logistic) (take 25))]
  (doseq [rf ["RF00114" "RF00521" "RF00555"]]
    (time (doall (build-sccs-pos-neg-sets
                  (->> (fs/directory-files
                        (fs/join "/home/kaila/Bio/Test/ROC-data/Sim" rf "CSV" )
                        ".cmsearch.csv")
                       first))))))

;;; All Single Context Sim
(binding [sccs-cut-points (->> (range -0.2 6 0.2) (map logistic) (take 25))]
  (let [base "/home/kaila/Bio/Test/ROC-data/Sim"]
    (doseq [rf (sort (set/difference
                      (->> "/home/kaila/Bio/Test/ROC-data/Sim/RF*" fs/glob
                           (map fs/basename) set)
                      (set ["RF00114" "RF00521" "RF00555"])))]
      (when (fs/directory? (fs/join base rf "SCCS-ROC"))
        (fs/rename (fs/join base rf "SCCS-ROC")
                   (fs/join base rf "ORIG-SCCS-ROC")))
      (time (doall (build-sccs-pos-neg-sets
                    (->> (fs/directory-files
                          (fs/join "/home/kaila/Bio/Test/ROC-data/Sim"
                                   rf "CSV" ) ".cmsearch.csv")
                         first)))))))

;;; Old, New
(let [all (->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
               (map #(->> % fs/basename (str/split #"-") first)))
      new (->> (fs/glob "/data2/BioData/RFAM-NC/NC/*.sto")
               (map #(->> % fs/basename (str/split #"-") first)))
      old (set/difference (set all) (set new))]
  (sort old))


(->> (fs/directory-files
      "/home/kaila/Bio/Test/ROC-data/Sim/RF00521/CLU-G" ".sto")
     (map #(sccs/save-hit-context-delta %)))


(defn sccs-roc-data [& {:keys [base RNAs]
                        :or {base "/home/kaila/Bio/Test/ROC-data"
                             RNAs ["L20" "S15" "S4" "S6S" "S6M"]}}]
  (for [rna RNAs]
    (let [dir (fs/join base rna "SCCS-ROC")
          cutpts (->> (fs/glob (fs/join dir "/sccs*-pos.ent"))
                      (map fs/basename)
                      (map #(->> % (str/split #"-") second)) sort)
          totdir (fs/join base rna "TotPosNeg")
          tot-pos (fs/join totdir (str rna "-pos.ent"))
          tot-neg (fs/join totdir (str rna "-neg.ent"))
          P (-> tot-pos get-entries count)
          N (-> tot-neg get-entries count)]
      (for [cpt cutpts
            :let [pos-ent (fs/join dir (str "sccs-" cpt "-pos.ent"))
                  neg-ent (fs/join dir (str "sccs-" cpt "-neg.ent"))]]
        (let [tp (count (tools/entry-file-intersect true pos-ent tot-pos))
              fn (count (tools/entry-file-intersect true neg-ent tot-pos))
              fp (count (tools/entry-file-intersect true pos-ent tot-neg))
              tn (count (tools/entry-file-intersect true neg-ent tot-neg))]
          [rna cpt
           {:TPR (true-positive-rate tp fn)
            :FPR (false-positive-rate fp tn)
            :MCC (mcc tp tn fp fn)
            :ACC (acc tp tn P N)
            :ppv (ppv tp fp)}])) )))

(def sccs-ecoli
     (sccs-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Ecoli"
      :RNAs ["S1" "S4" "S7" "S8" "S15"]))


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
            plot-dir (fs/join base "Plots")
            chart-file (fs/join plot-dir
                                (str (nth RNAs i) "-roc-plot.png"))]
        (when (not (fs/directory? plot-dir))
          (fs/mkdir plot-dir))
        (incanter.core/view chart)
        (incanter.core/save
         chart chart-file :width 500 :height 500)))))

(build-roc-curves)


(defn csv-roc-data
  [& {:keys [base RNAs subname]
      :or {base "/home/kaila/Bio/Test/ROC-data"
           RNAs ["L20" "S15" "S4" "S6S" "S6M"]
           subname ""}}]
  (let [subname (if (= subname "") "-" (str "-" subname "-"))
        sccs-data (sccs-roc-data :base base :RNAs RNAs)
        eval-data (eval-roc-data :base base :RNAs RNAs)
        step (/ 1.0 (->> sccs-data first count))
        cols ["RNA" "Cutpt" "TPR" "FPR" "Type"]
        guess (map #(do ["" "" (str %) (str %) "Guess"]) (range 0.0 1.0 step))
        xf #(map (fn[[nm cpt vm]]
                   (let [cpt (if (= %1 "SCCS")
                               (/ (floor (* (Double. cpt) 10000.0)) 10000.0)
                               cpt)]
                   (cons nm (map str [cpt (vm :TPR) (vm :FPR) %1]))))
                 %2)
        data-dir (fs/join base "ROCdata")]
    (when (not (fs/directory? data-dir)) (fs/mkdirs data-dir))
    (println step guess)
    (doseq [i (range 0 (count RNAs))]
      (let [sccsi (xf "SCCS" (nth sccs-data i))
            evali (xf "EVAL" (nth eval-data i))
            all (concat sccsi evali guess)]
        (spit (fs/join data-dir (str (ffirst all) subname "both.csv"))
              (csv/write-csv (cons cols all)))))))

(csv-roc-data :base "/home/kaila/Bio/Test/ROC-data/Ecoli"
              :RNAs ["L1" "L10" "L20" "L4" "S1" "S15" "S2" "S4" "S7" "S8"])





(def bioinfo-table-template
     ["{\\scriptsize \\begin{tabular}{ll|lllll}\\toprule"
      "ncRNA  & xxctpt  & AUC    & TPR    & SPC    & PPV    & MCC " "\\toprule"
      "S1EC   & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "S4EC   & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "S7EC   & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "S8EC   & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "S15EC  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "L20BS  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "S4BS   & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "S6BS   & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "S15BS  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam1  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam2  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam3  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam4  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam5  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam6  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam7  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam8  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam9  & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam10 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam11 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam12 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam13 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam14 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam15 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam16 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam17 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam18 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam19 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "Rfam20 & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "\\botrule"
      "\\end{tabular}}"])



(defn build-sim-config-files
  []
  (let [sim-dir "/home/kaila/Bio/Test/ROC-data/Sim"
        dirs (->> (fs/glob (fs/join sim-dir "RF0*"))
                  sort)
        runtxts (->> (fs/directory-files sim-dir ".txt")
                     (map #(do [(fs/basename %) (io/read-lines %)])))]
    (doseq [d dirs]
      (let [dnm (fs/basename d)
            xpart (str/lower-case dnm)
            sto (->> d (#(fs/directory-files % ".sto"))
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

(build-sim-config-files)






;;;


#_(defn get-rf-full-sets-cmcsvs [& {:keys [cutpt] :or {cutpt 1.0e-12}}]
  (->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
       (map #(-> % (fs/join "CSV/*cmsearch.csv") fs/glob first))
       (map #(do [(->> % fs/basename (str/split #"-") first)
                  (-> % slurp csv/parse-csv  rest butlast)]))
       (sort-by first)
       (map (fn[[n recs]] [n (sort-by (fn[x] (Float. (nth x 9))) recs)]))
       (map (fn[[n recs]] [n (take-until #(> (Float. (nth % 9)) cutpt) recs)]))
       #_(map (fn[[n recs]] [n (count recs)]))))



(defn csv-recs-to-entries [csv-info-recs]
  (keep (fn[[nm s e sd]]
          (when (fs/exists? (fs/join (default-genome-fasta-dir nm)
                                     (str nm ".fna")))
            (str nm "/"
                 (if (> (Integer. s) (Integer. e))
                   (str e "-" s "/-1")
                   (str s "-" e "/1")))))
        csv-info-recs))

(defn sane-csv-info [rows]
  (reduce (fn[res rec]
            (conj res [(rec 0) (rec 3) (rec 4) (rec 9)]))
          [] rows))

(defn get-rf-rs58-sets [& {:keys [cutpt] :or {cutpt 1.0e-1}}]
  (->> (fs/glob "/home/kaila/Bio/Test/ROC-data/TMP/Sim/RF*")
       (map #(-> % (fs/join "Pre*") fs/glob first))
       (map #(-> % (fs/join "*.csv") fs/glob first))
       (map #(do [(->> % fs/basename (str/split #"-") first)
                  (-> % slurp csv/parse-csv rest butlast sane-csv-info)]))
       (sort-by first)
       (map (fn[[n recs]] [n (sort-by (fn[x] (Float. (last x))) recs)]))
       (map (fn[[n recs]] [n (take-until #(> (Float. (last %)) cutpt) recs)]))
       (map (fn[[n recs]] [n (->> recs csv-recs-to-entries
                                 gen-name-seq-pairs vec)]))
       #_(map (fn[[n recs]] [n (count recs)]))))


(defn get-rf-full-sets []
  (->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
       (map #(-> % fs/basename (str "-full.sto")))
       (map #(fs/join "/data2/Bio/RFAM" %))
       (map #(do [(->> % fs/basename (str/split #"-") first)
                  (-> % (read-seqs :info :both))]))
       (sort-by first)
       (map (fn[[n sqs]]
              [n (vec (map (fn[[ent sq]]
                             [ent (->> sq degap-seqs norm-elements)])
                           sqs))]))
       #_(map (fn[[n recs]] [n (count recs)]))))


(defn get-set-counts [& {:keys [cpt58]
                         :or {cpt58 1.0e-1}}]
  (let [rs58s (get-rf-rs58-sets :cutpt cpt58)
        rfFulls (get-rf-full-sets)]
    (map (fn[[n58 recs58] [nfull recsfull]]
           (let [cnt58 (count recs58)
                 cntfull (count recsfull)
                 m (min cnt58 cntfull)
                 mx (max cnt58 cntfull)
                 cnt (if (not= 0 m) m (min mx 100))
                 cnt (cond (< cnt 20) (if (< mx 100) mx 100)
                           (> cnt 100) 100
                           :else cnt)]
             {:name n58 :cnt cnt :cnt2 (* 2 cnt)
              :cnt58 cnt58 :cntfull cntfull}))
         rs58s rfFulls)))

(defn get-sample-sets [& {:keys [cpt58]
                          :or {cpt58 1.0e-1}}]
  (let [rs58s (get-rf-rs58-sets :cutpt cpt58)
        rfFulls (get-rf-full-sets)]
    (map (fn[[n58 recs58] [nfull recsfull]]
           (let [cnt58 (count recs58)
                 cntfull (count recsfull)
                 m (min cnt58 cntfull)
                 mx (max cnt58 cntfull)
                 cnt (if (not= 0 m) m (min mx 100))
                 cnt (cond (< cnt 20) (if (< mx 100) mx 100)
                           (> cnt 100) 100
                           :else cnt)]
             {:name n58 :cnt cnt :cnt58 cnt58 :cntfull cntfull
              :set58 (take cnt recs58)
              :setfull (take (* 2 cnt) recsfull)}))
         rs58s rfFulls)))


(def sim-sample-sets (get-sample-sets))

(->> sim-sample-sets (map :set58)
     (map #(do [(count %) (mean (map (fn[[ent sq]] (count sq)) %))])))
(->> sim-sample-sets (map :setfull)
     (map #(do [(count %) (mean (map (fn[[ent sq]] (count sq)) %))])))

(let [base "/home/kaila/Bio/Test/ROC-data/Sim/AllDBSeeds"
      samples-58 (->> sim-sample-sets (map #(do [(:name %) (:set58 %)])))
      samples-full (->> sim-sample-sets (map #(do [(:name %) (:setfull %)])))
      sample-sets (map (fn[[n58 ents-sqs58] [nfull ents-sqsfull]]
                         (if (not= n58 nfull)
                           (raise :name-mismatch :n58 n58 :nfull nfull)
                           [n58 (concat ents-sqs58 ents-sqsfull)]))
                       samples-58 samples-full)]
  (doseq [[nm ents-sqs] sample-sets]
    (nms-sqs->fasta-file ents-sqs (fs/join base (str nm "-sample.fna")))))



;;; Generate simulated context locations.
;;;
(def exon-select
     "select sfqv.term_id,sfqv.value,loc.start_pos,loc.end_pos,loc.strand
             from bioentry as be,
                  seqfeature as sf,
                  seqfeature_qualifier_value as sfqv,
                  location as loc
             where be.bioentry_id=sf.bioentry_id and
                   sf.seqfeature_id=sfqv.seqfeature_id and
                   sf.type_term_id=12 and sfqv.term_id=14 and
                   loc.seqfeature_id=sf.seqfeature_id and
                   be.name=\"/name/\" and
                   loc.start_pos > /start/ and
                   loc.strand = /strand/ limit 30")
"NC_017986"

(->>
 (db/sql-query
  (->> exon-select
       (str/replace-re #"/name/" "NC_018643")
       (str/replace-re #"/start/" "1341411")
       (str/replace-re #"/strand/" "1")))
 ensure-vec (drop 1)
 (map #(% :start_pos)))


(defn get-sim-entry [entry]
  (let [[n [s e] st] (entry-parts entry)
        sz (abs (- e s))
        starts (->> (db/sql-query
                     (->> exon-select
                          (str/replace-re #"/name/" n)
                          (str/replace-re #"/start/" (str s))
                          (str/replace-re #"/strand/" (str st))))
                    ensure-vec (drop 1)
                    (map #(% :start_pos)))]
    (set (map #(let [s' (- % sz 50)
                     e' (+ s' sz)]
                 (make-entry n s' e' st))
              starts))))

(defn gen-sim-entries [in-files]
  (let [ot-files (->> in-files (map #(str/replace-re #"seed" "simout" %))
                      (map #(fs/replace-type % ".ent")))
        entry-sets (->> in-files
                        (map #(read-seqs % :info :name))
                        (map #(apply set/union (map get-sim-entry %))))]
    (doseq [[f entries] (partition-all 2 (interleave ot-files entry-sets))]
      (gen-entry-file entries f))))

(gen-sim-entries (fs/glob "/data2/BioData/RFAM-NC/NC/*seed*.sto"))


(->> (get-csv-entry-info "/home/kaila/Bio/Test/ROC-data/Sim/RF00162/PrelimCSV/RF00162-seed-NC.sto.aggregate.cmsearch.csv")
     (filter #(> (Double. (nth % 3)) 1.0e-2))
     (map (fn[[n s e ev]]
            (let [[s e] [(Integer. s) (Integer. e)]
                  [s e std] (if (> s e) [e s -1] [s e 1])]
              [(make-entry n s e std) (Double. ev)])))
     (sort-by second)
     (take 10))
(map #(->> % gen-name-seq-pairs (map second)))

(defn gen-sim-db []
  (let [sim-files (sort (fs/glob "/data2/BioData/RFAM-NC/NC/*simout*.ent"))
        simctxs (->> sim-files
                     (map #(read-seqs % :info :name))
                     (map (fn[ents]
                            (map #(let [[nm [s e] std] (entry-parts %)
                                        rdelta (+ (- e s) 850)]
                                    (gen-name-seq % :ldelta 100 :rdelta rdelta))
                                 ents))))
        fullsamps (->> "/home/kaila/Bio/Test/ROC-data/Sim/AllDBSeeds/*.fna"
                       fs/glob sort
                       (map #(vec (read-seqs % :info :data))))]

    (println :sim-files (count sim-files)
             :simctxs (count simctxs)
             :fullsamps (count fullsamps))

    (doseq [[f nm-sqs] (map (fn[sim-file seed-simctxs ncRNAs]
                              [(fs/replace-type sim-file ".fna")
                               (map (fn[[ent ctxseq]]
                                      (let [ncRNA (rand-nth ncRNAs)
                                            [nm [s e] std] (entry-parts ent)
                                            orig-rna-gap-sz (- (- e s) 950)
                                            drop-amt (+ 100 orig-rna-gap-sz)
                                            suffix (str/drop drop-amt ctxseq)
                                            new-ent (make-entry
                                                     nm s
                                                     (dec (+ s 100 (count ncRNA)
                                                             (count suffix)))
                                                     std)]
                                        [new-ent
                                         (str (str/take 100 ctxseq)
                                              ncRNA
                                              (str/drop drop-amt ctxseq))]))
                                    seed-simctxs)])
                            sim-files (doall simctxs) (doall fullsamps))]
      (io/with-out-writer f
        (doseq [[ent sq] nm-sqs]
          (println (str ">" ent))
          (println sq))))))


(defn gen-sim-db-0 []
  (let [sim-files (sort (fs/glob "/data2/BioData/RFAM-NC/NC/*simout*.ent"))
        simctxs (->> sim-files
                     (map #(read-seqs % :info :name))
                     (map (fn[ents]
                            (map #(let [[nm [s e] std] (entry-parts %)
                                        rdelta (+ (- e s) 850)]
                                    (gen-name-seq % :ldelta 100 :rdelta rdelta))
                                 ents))))
        fullsamps (->> sim-files (map fs/basename)
                       (map #(str/replace-re #"simout-NC.ent$" "full.sto" %))
                       (map #(fs/join "/data2/Bio/RFAM" %))
                       (map #(read-seqs % ))
                       (map degap-seqs) (map norm-elements)
                       (map set) (map vec))]
    (doseq [[f nm-sqs] (map (fn[sim-file seed-simctxs ncRNAs]
                              [(fs/replace-type sim-file ".fna")
                               (map (fn[[ent ctxseq] ncRNA]
                                      (let [[nm [s e] std] (entry-parts ent)
                                            orig-rna-gap-sz (- (- e s) 950)
                                            drop-amt (+ 100 orig-rna-gap-sz)
                                            suffix (str/drop drop-amt ctxseq)
                                            new-ent (make-entry
                                                     nm s
                                                     (dec (+ s 100 (count ncRNA)
                                                             (count suffix)))
                                                     std)]
                                        [new-ent
                                         (str (str/take 100 ctxseq)
                                              ncRNA
                                              (str/drop drop-amt ctxseq))]))
                                    seed-simctxs
                                    ncRNAs)])
                            sim-files (doall simctxs) (doall fullsamps))]
      (io/with-out-writer f
        (doseq [[ent sq] nm-sqs]
          (println (str ">" ent))
          (println sq))))))








;;; Backup original search good output CSVs
;;;
(doseq [f (->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
               (map #(fs/join % "CSV"))
               (map #(first (fs/directory-files % ".cmsearch.csv"))))]
  (fs/copy f (fs/replace-type f ".csv.bkup")))


;;; Union the original seeds (in csv format - see below) onto original
;;; search output CSVs
;;;
(let [dirs (->> (fs/directory-files "/home/kaila/Bio/Test/ROC-data/Sim" "")
                (filter #(re-find #"RF0" %)))
      csvs (->> dirs (map #(-> % fs/basename (str "-seed-NC.csv")))
                 (map #(fs/join %1 "CSV" %2) dirs))
      CSVs (map #(fs/join % "CSV") dirs)
      out-csvs (map #(first (fs/directory-files % "cmsearch.csv.bkup")) CSVs)
      csv-pairs (sort-by first (partition-all 2 (interleave out-csvs csvs)))]
  csv-pairs)
  (doseq [[ocsv csv] csv-pairs]
    (io/with-out-writer (fs/replace-type ocsv "")
      (doseq [l (io/read-lines ocsv)] (println l))
      (doseq [l (io/read-lines csv)] (println l)))))


;;; Create seed CSVs to be unioned to output CSV results
;;;
(defn gen-sto-csvs [& {:keys [base RNAs name]
                       :or {base "/home/kaila/Bio/Test/ROC-data/Ecoli"
                            RNAs ["S1" "S2" "S4" "L4" "L1" "L10"]
                            name "foo.sto"}}]
  (let [seeds (map #(fs/join base % name) RNAs)
        csvs (map #(-> base (fs/join % (str name "-missing.csv"))) RNAs)]
    (doseq [[seed csv] (partition-all 2 (interleave seeds csvs))]
      (->> (get-sto-as-csv-info seed :ev 1.3)
           (map (fn[[nm s e ev _ _ std]]
                  (csv/csv-to-stg
                   (map str [nm 1 2 s e std 11 22 0.0 ev 0.0 0.0]))))
           doall (str/join "\n")
           (spit csv)))))



#=GF CTXSZ 260

(->> (get-csv-entry-info "/home/kaila/Bio/Test/ROC-data/Sim/RF00162/CSV/RF00162-seed-NC.sto.all-sim-db.cmsearch.csv")
     (map (fn[[n s e ev]]
            (let [[s e] [(Integer. s) (Integer. e)]
                  [s e std] (if (> s e) [e s -1] [s e 1])]
              [(make-entry n s e std) (Double. ev)])))
     (sort-by second)
     (take 10))



(->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
     #_(filter #(not (fs/exists? (fs/join % "SCCS-ROC"))))
     sort (map #(let [x (fs/join % "CLU") y (fs/join % "CLU-G")]
                  (if (fs/directory? x) x y)))
     (map #(do [(->> % fs/split butlast last)
                (slurp (fs/join % "delta.txt"))
                (fs/directory-files % ".sto")]))
     (map #(do [(first %)
                (second %)
                (map (fn[f] (count (read-seqs f :info :name))) (last %))]))
     (map (fn [[n d x]] [n d (sum x) x])))

(->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
     (filter #(not (fs/exists? (fs/join % "SCCS-ROC"))))
     sort (map #(fs/join % "CLU"))
     (map #(do [(->> % fs/split butlast last) (fs/directory-files % ".sto")]))
     (map #(do [(first %)
                (map (fn[f] (count (read-seqs f :info :name))) (second %))]))
     (map (fn [[n x]] [n (sum x) x])))

(->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
     sort (map #(fs/join % "CLU"))
     (map #(fs/directory-files % ".sto"))
     (map #(map (fn[f] (count (read-seqs f :info :name))) %))
     (map #(do [(sum %) %]))
     (interleave (->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
                      sort
                      (map fs/basename)))
     (partition-all 2))

(->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
     sort (map #(first (fs/directory-files % ".sto")))
     (map #(do [(->> % fs/dirname fs/basename)
                (-> % fs/dirname (fs/join "CLU") (fs/join "delta.txt") slurp)
                (count (read-seqs % :info :name))]))
     (sort-by #(->> % second (str/drop 1) Integer.)))


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