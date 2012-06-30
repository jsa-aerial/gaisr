;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                         P O S T - D B - C S V                            ;;
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

(ns edu.bc.bio.gaisr.post-db-csv

  "Post processing of pipeline results.  Revolves around filtering
   output (remove duplicates, ev cutoffs, etc) as well as formatting
   to DB table (and corresponding CSV) formats"

  (:require [clojure.contrib.sql :as sql]
            [org.bituf.clj-dbcp :as dbcp]
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure-csv.core :as csv]
            [clojure.contrib.json :as json]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.xml :as xml]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        [edu.bc.log4clj
         :only [create-loggers log>]]

        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.tools

        [edu.bc.bio.gaisr.db-actions
         :only [+start-delta+ base-info-query hit-features-query]]

        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)])

  (:import javax.sql.DataSource
           com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource
           [java.util StringTokenizer]
           [java.io StreamTokenizer]))


(def test-file
     "/data2/Bio/Work/ECS15_032811-2009-12-17-utrs-ribs-refseq32-microbial.cmzasha.csv")


(def +base-codepoint+ (first (str/codepoints "a")))

(def col-index-map
     (reduce (fn[m [l cp]] (assoc m (keyword (str l)) (- cp +base-codepoint+)))
             {}
             (map #(do [%1 %2])
                  "abcdefghijklmnopqrstuvwxyz"
                  (str/codepoints "abcdefghijklmnopqrstuvwxyz"))))




;;; GAISR CSV Header
;;;
(def +cmsearch-csv-header+
     "gaisr name,orig-start,orig-end,hit-start,hit-end,hit-strand,hit-rel-start,hit-rel-end,score,evalue,pvalue,gc,structure,tgt-seq,orig-tgt-seq")


(defn legacy-csv [rows]
  (let [newnes [:rfam :overlap :new]]
    (loop [entries []
           rows (drop 1 rows)]
      (let [entry (first rows)]
        (if (< (count entry) 24)
          (sort #(string-less? (%1 0) (%2 0)) entries)
          (recur (conj entries
                       [(entry 4) (entry 9) (entry 10) (entry 24)
                        (entry 11)
                        (newnes (Integer/parseInt (entry 15)))
                        (entry 13)])
                 (drop 1 rows)))))))

(defn gaisr-csv [rows]
  (loop [entries []
         rows (drop 1 rows)]
    (let [entry (first rows)]
      (if (< (count entry) 14) ; bogus due to csv/parse-csv....
        (sort #(string-less? (%1 0) (%2 0)) entries)
        (recur (conj entries
                     [(entry 0) (entry 3) (entry 4) (entry 9)
                      (entry 8)
                      :new "N/A"])
               (drop 1 rows))))))

(defn get-entries [csv-hit-file]
  (let [file (fs/fullpath csv-hit-file)
        rows (csv/parse-csv (slurp file))
        head (first rows)]
    ;;(prn head)
    ;;(prn (head 4) (head 9) (head 10) (head 24) (head 11) (head 15) (head 13))
    (if (= (first head) "gaisr name")
      (gaisr-csv rows)
      (legacy-csv rows))))


(defn filter-entries [entries]
  (keep #(when (and (= (% 5) :new)
                    (not (re-find #"^NZ_" (% 0)))
                    (not (re-find #"duplicate" (% 6))))
           %)
        entries))


(defn efile->names [filespec]
  (with-local-vars [names #{}]
    (do-text-file [(fs/fullpath filespec)]
      (when (and (> (count $line) 3)
                 (or (= ">gi"  (subs $line 0 3))
                     (= "#=GS" (subs $line 0 4))
                     (= "NC_"  (subs $line 0 3))))
        (var-set names (conj (var-get names)
                             (first (re-find #"N(C|S|Z)_[0-9A-Z]+" $line))))))
    (var-get names)))

;;; ---------------- Collect matching and overlapping N* -------------------

(defn make-start-set [filespec]
  (filter #(or (= :new (% 1))
               (= :rfam (% 1))
               (= :overlap (% 1))
               (re-find #"duplicate" (% 2)))
          (map #(let [nm (first (str/split #"\." (% 0)))]
                  [nm (% 5) (or (first (re-find #"N(C|Z|S)_[0-9A-Z]+" (% 6)))
                                nm)])
               (get-entries filespec))))

(defn make-qset [start-set]
  (reduce (fn[m [n* _ mn*]]
            (if (= (subs mn* 0 1) "N")
              (assoc m mn* (conj (get m mn* #{}) n*))
              m))
          {} start-set))

(defn make-entry-hit-map [start-set]
  (reduce (fn[m [n* _ mn*]]
            (assoc m
              n* (conj (get m n* #{})
                       (if (= (subs mn* 0 1) "N") mn* :dvo))))
          {} start-set))

(defn entry->eqv-class [entry hit-map qset]
  (let [hits (hit-map entry)]
    (reduce (fn[s h] (set/union s (if (not= h :dvo) (qset h) #{})))
            #{} hits)))

(defn efile-csvhits->names-and-matches [efile csv-hit-file & [out-filespec]]
  (let [start-set (make-start-set csv-hit-file)
        qset (make-qset start-set)
        hit-map (make-entry-hit-map start-set)
        names (efile->names efile)
        name-ecs (map #(do [% (entry->eqv-class % hit-map qset)]) names)
        all-names (reduce (fn[s [n ec]] (set/union s ec))
                          #{} name-ecs)]
    (when out-filespec
      (io/with-out-writer (fs/fullpath out-filespec)
        (doseq [[n ecs] name-ecs]
          (println n "-->" ecs))
        (println "\n****\n")
        (doseq [n all-names]
          (println n))))
    (list name-ecs all-names)))


;;; ----------------- End matching and overlapping N* ----------------------


;;; START FUGLY HACK---------------------------------------------------------
;;;
(defn twiddle-name
  "If name = prev-name, create a 'versioned' name by appending /n to
   name.  n is an incrementing count that starts at 2.  So, generates
   a sequenceof names: name, name/2, name/3, ...

   This nonsense has to be accounted for when generating entry files
   from this information in actions/action-request :genfasta and kin.
  "
  [name prev-name]
  (if (= name prev-name)
    (let [i (re-find #"/[0-9]" name)
          i (when i (Integer. (subs i 1)))]
      (if (not i)
        (str name "/2")
        (str (str/take (dec (count name)) name) (inc i))))
    name))

(defn base-info
  "A hack for ensuring proper alignment of hit feature info, hit loci
   and the base info for a hit when there are multiple hits in a
   species or the database is missing the species (due to unloaded
   data source or some other anomaly)"
  [names]
  (loop [prev-name (first names)
         ns names
         infos (base-info-query names)
         ninfos []]
    (if (empty? ns)
      ninfos
      (let [n (first ns)
            n (subs n 1 (dec (count n)))
            info-name (:name (first infos))
            gbid      (:gbid (first infos))]
        (if (not= n info-name)
          (recur (twiddle-name n prev-name)
                 (drop 1 ns)
                 infos
                 (conj ninfos {:name (twiddle-name n prev-name)
                               :version 0 :description "NA"
                               :bioentry_id (gen-kwuid) :gbid gbid
                               :taxon_id 0 :taxname "NA" :ancestors "NA"
                               :delta +start-delta+}))
          (recur n
                 (drop 1 ns)
                 (drop 1 infos)
                 (conj ninfos (first infos))))))))
;;;
;;; END FUGLY HACK-----------------------------------------------------------

(defn gather-hit-features [hits]
  (let [hfs (map (fn[[name start end & tail]]
                   (let [name (first (str/split #"\." name))
                         start (Integer. start)
                         end (Integer. end)
                         tmp start
                         start (if (< end start) end start)
                         end (if (= start end) tmp end)]
                     (hit-features-query name start end)))
                 hits)
        bi (base-info
            (map #(cl-format nil "~S" (first (str/split #"\." (first %1))))
                 hits))
        hit-loci-ev (map (fn[[_ start end ev & tail]]
                           (let [start (Integer. start)
                                 end (Integer. end)
                                 tmp start
                                 strand (if (< end start) -1 1)
                                 start (if (< end start) end start)
                                 end (if (= start end) tmp end)]
                             [start end ev strand]))
                         hits)]
    (map (fn [bimap [s e ev st] hf]
           (let [hit-pseudo-feature
                 {:sfid 0 :sftype "hit"
                  :locs (seq [{:start s :end e :len (- e (dec s)) :strand st}])
                  :nvs ()}]
           (assoc bimap
             :hitloc (str s "-" e)
             :hitstrand st
             :evalue ev
             :hit_features (conj hf hit-pseudo-feature)
             :sfcount (count hf)
             :delta +start-delta+)))
         bi hit-loci-ev hfs)))


(defn process-hit-file [filespec]
  (let [f (fs/fullpath filespec)
        hitseq (-> f get-entries filter-entries)]
    (if (empty? hitseq)
      {:error (cl-format nil "File '~A' has no nonredundant content" filespec)}
      (gather-hit-features hitseq))))




;;; ----- Count various frequencies of characteristics of output (from csvs)


(defn freq [csv-dir cnt-fn & args]
  (let [files (sort (fs/directory-files csv-dir ".cmsearch.csv"))
        names (map #(second (str/split #"\." (re-find #"sto\..*\.cmsearch" %)))
                   files)]
    (map #(do [%1 (apply cnt-fn %2 args)])
         names files)))


(defn ev-freq [cmsearch-csv ev-cutoff]
  (reduce (fn [n x]
            (let [ev (Float. (x 9))]
              (if (< ev ev-cutoff) (inc n) n)))
          0 (butlast (drop 1 (csv/parse-csv
                              (slurp (fs/fullpath cmsearch-csv)))))))

(defn ev-freq-ss [dirdir ev-cutoff outss-filespec]
  (let [ev-cnts (fs/dodir
                 dirdir
                 #(butlast (drop 1 (sort (filter fs/directory?
                                                 (fs/directory-files % "")))))
                 #(do [(fs/basename %) (freq % ev-freq ev-cutoff)]))
        cols (csv/csv-to-stg (cons "Names" (map first ev-cnts)))
        rows (apply map vector (cons (map first (second (first ev-cnts)))
                                     (map #(map second (second %)) ev-cnts)))
        rows (csv/write-csv (map #(map str %) rows))]
    (io/with-out-writer (fs/fullpath outss-filespec)
      (println cols)
      (print rows))))




;;; ----- Generate positive/negative kernel training sets from CSVs -----


(defn get-positive-negative-sets
  "Generate positive and negative training sets from cmsearch out
  CSVs.  See ...gaisr.pipeline/cmsearch-out-csv, et.al.  Take the
  found hits and separate them by EV-CUTOFF.  EV-CUTOFF is either a
  single floating point value or a two element vector of floating
  point values denoting a range.  It defaults to 1.0.

  Hits ->
    Positives {h | (< (evalue h) ev-cutoff), h in Hits}
    Negatives {h | (<= ev-cutoff (evalue h)), h in Hits}

  Hits, ev-cutoff = [s e]->
    Positives {h | (< (evalue h) s), h in Hits}
    Negatives {h | (< e (evalue h)), h in Hits}"

  [cmsearch-out-csv & {ev-cutoff :ev-cutoff :or {ev-cutoff 1.0}}]

  (let [[hd & rows] (csv/parse-csv (slurp cmsearch-out-csv))
        evs (if (coll? ev-cutoff) (first ev-cutoff) ev-cutoff)
        eve (if (coll? ev-cutoff) (second ev-cutoff) (float Integer/MAX_VALUE))
        nmpos 0
        epos (first (seq/positions #(= "evalue" %) hd))
        seqpos (first (seq/positions #(= "orig-tgt-seq" %) hd))
        [p1 n1] (seq/separate #(and (> (count %) 1) ; Bogus empty set check
                                    (< (Float. (nth % epos)) evs))
                              rows)
        n1 (filter #(and (> (count %) 1) ; Bogus empty set check
                         (< eve (Float. (nth % epos))))
                   n1)
        p1 (map #(do [(nth % nmpos) (nth % seqpos)]) p1)
        n1 (map #(do [(nth % nmpos) (nth % seqpos)]) n1)]
    [p1 n1]))


(defn gen-positive-negative-sets
  "Generate positive and negative training sets for various kernel/SVM
   from previously proposed candidates (via cmsearch).  Take the
   previous candidates in CSV format (see
   ...gaisr.pipeline/cmsearch-out-csv, et.al) and separate them by
   evalue.

   EV-CUTOFF is either a single floating point value (such as 0.1) or
   a two element vector of floating point values denoting the range of
   evalues for negative examples.  Defaults to single value 1.0.

   If EV-CUTOFF single value those entries with evalues < ev-cutoff
   are placed, in fasta format, in POS-OUT-FSPEC file and those with
   evalues >= ev-cutoff are placed, in fasta format, in NEG-OUT-FSPEC.

   If EV-CUTOFF is a range [s e] then those entries with evalues < s
   are placed, in fasta format, in pos-out-fspec and those enties with
   s <= evalue < e are placed, in fasta format, in neg-out-fspec."

  [cmsearch-out-csv pos-out-fspec neg-out-fspec
   & {ev-cutoff :ev-cutoff maxseqs :maxseqs
      :or {ev-cutoff 1.0 maxseqs 3319}}]

  (let [[pset nset] (get-positive-negative-sets
                     cmsearch-out-csv :ev-cutoff ev-cutoff)
        pset (if (> (count pset) maxseqs) (take maxseqs (shuffle pset)) pset)
        nset (if (> (count nset) maxseqs) (take maxseqs (shuffle nset)) nset)
        pout (fs/fullpath pos-out-fspec)
        nout (fs/fullpath neg-out-fspec)]
    (nms-sqs->fasta-file pset pout)
    (nms-sqs->fasta-file nset nout)
    [pout nout]))


(defn get-combined-csv
  "Build a file of combined csv content.  CSV-OR-CSV-DIR is a filespec
  of either a single csv file or a directory containing csv files,
  either in subdirs or directly.

  PRED is a predicate filter of a single filespec argument. Only those
  files passing it are returned.  If given, PRED should return its
  input argument, or nil (it acts as a filter to KEEP).

  If CSV-OR-CSV-DIR is a single csv file, return it if pred returns
  true on it or nil otherwise.

  If CSV-OR-CSV-DIR is a collection, run filter on it with pred and
  return a temp file containing the combined entries for all elements
  that pass pred.

  If it is a directory, then if SUBDIRS is true (default) the csvs are
  in subdirs of the given directory, otherwise they are taken to be
  directly in the dir.  Return a temp file of the combined entries for
  all such files which pass pred."

  [csv-or-csv-dir
   & {pred :pred subdirs :subdirs :or {pred identity subdirs true}}]

  (let [f-or-d (if (coll? csv-or-csv-dir)
                 (map fs/fullpath csv-or-csv-dir)
                 (fs/fullpath csv-or-csv-dir))
        reducer (fn [csvs]
                  (reduce
                   (fn[om f]
                     (reduce
                      (fn[m l]
                        (let [k (ffirst (csv/parse-csv l))]
                          (assoc m k l)))
                      om
                      (drop 1 (io/read-lines f))))
                   {} csvs))]
    (cond
     (coll? f-or-d)
     (let [tmp-file (fs/tempfile "catcsvs-" ".csv")]
       (io/with-out-writer tmp-file
         (println +cmsearch-csv-header+)
         (let [unique-line-map (reducer (filter pred f-or-d))]
           (doseq [l (vals unique-line-map)]
             (println l))))
       tmp-file)

     (not (fs/directory? f-or-d))
     f-or-d ; file - assumes it is a csv

     :Else ;  directory
     (let [tmp-file (fs/tempfile "catcsvs-" ".csv")]
       (io/with-out-writer tmp-file
         (println +cmsearch-csv-header+)
         (let [base f-or-d
               dirs (if subdirs ; csvs in subdirs of given dir
                      (keep #(let [f (fs/join base %)]
                               (when (fs/directory? f) f))
                            (fs/listdir base))
                      ;; else csvs in given dir
                      [base])]
           (doseq [d dirs]
             (let [csvs (keep #(when (and (re-find #"cmsearch.csv$" %)
                                          (pred %))
                                 (fs/join d %))
                              (fs/listdir d))
                   unique-line-map (reducer csvs)]
               (doseq [l (vals unique-line-map)]
                 (println l))))))
       tmp-file))))


(defn gen-aligned-pos-neg-sets

  [cm csv-or-csv-dir base-out
   & {pred :pred subdirs :subdirs ev-cutoff :ev-cutoff
      :or {pred identity subdirs false ev-cutoff 0.00001}}]

  (let [cm (fs/fullpath cm)
        bos (map #(fs/fullpath (str base-out % ".x"))
                 ["-pos" "-neg"])
        [pfna nfna] (map #(fs/replace-type % ".fna") bos)
        [psto nsto] (map #(fs/replace-type % ".sto") bos)]

    (gen-positive-negative-sets
     (get-combined-csv csv-or-csv-dir :subdirs subdirs :pred pred)
     pfna nfna :ev-cutoff ev-cutoff)
    (map (fn [seqf stof]
           (if (fs/empty? seqf)
             (do (io/with-out-writer stof
                   (println (str "# GEN-ALIGNED-POS-NEG-SETS: SeqFile '"
                                 seqf "' is empty!")))
                 stof)
             (cmalign cm seqf stof)))
         [pfna nfna]
         [psto nsto])))



(defn gen-aligned-training-sets [dir & args]
  (let [cm (first (fs/directory-files dir ".cm"))
        csv-groups (vals (group-by
                          (fn[f]
                            (->> f (str/split #"sto\.") second
                                 (str/split #"-") first
                                 (str/split #"\.") last))
                          (fs/directory-files dir ".cmsearch.csv")))
        bouts (map #(str/replace-re #"[0-9]*$" ""
                                    (fs/replace-type (first %) ["" ""]))
                   csv-groups)]
    (map #(apply gen-aligned-pos-neg-sets cm %1 %2 args)
         csv-groups
         bouts)))


;;; (gen-aligned-training-sets
;;;  "/data2/Bio/Training/MoStos/ECS4" :ev-cutoff 0.00001)
