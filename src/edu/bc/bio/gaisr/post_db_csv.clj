;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                         P O S T - D B - C S V                            ;;
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
            [clojure.java.shell :as sh]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        [edu.bc.log4clj
	 :only [create-loggers log>]]
        edu.bc.bio.seq-utils
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
(defn twiddle-name [name prev-name]
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


