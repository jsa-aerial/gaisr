;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                          D B - D O W N L O A D S                         ;;
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

(ns edu.bc.bio.db-downloads

  "HTTP and FTP public DB crawlers and auto downloader.  Initially,
  setup and configured for NCBI RefSeq genome/fasta downloads and NCBI
  refseq/genome preformatted BLAST dbs"

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.contrib.properties :as prop]
            [clojure-csv.core :as csv]
            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.net-utils
        [edu.bc.bio.sequtils.files
         :only [sto-GC-and-seq-lines join-sto-fasta-file]]
        [edu.bc.log4clj
         :only [create-loggers log>]]

        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)])

  (:import [java.lang Thread]))


(create-loggers
 [[:logger "refseq-ftp"
   {:level :info :appenders "refftp"}]
  [:appender "refftp"
   {:type :rolling
    :filespec "./logs/ref-ftp.log"
    :max-version 10
    :max-size 1000000}]
  ])


(defn refseq-files [ftp-set]
  (keep #(let [name (.getName %)]
           (when (and (.isFile %)
                      (re-find #"[0-9]+\.genomic\.(gbff|fna)\.gz" name))
             name))
        ftp-set))


(defn nz-files [ftp-set]
  (keep #(let [name (.getName %)]
           (when (and (.isFile %)
                      (re-find #"microbialNZ.*\.genomic\.(gbff|fna)\.gz" name))
             name))
        ftp-set))


(defn list-refseq-files [url]
  (ftp-dir url refseq-files))


(defn retrieve-refseq-files [url files dir]
  (doseq [file files]
    (let [filespec (str dir "/" file)]
      (log> "refseq-ftp" :info
            "RefSeq ftp: Fetching ~A --> ~A" file filespec)
      (retrieve-file url file (str dir "/" file))
      (Thread/sleep 500))))


;;; "ftp://anonymous:user%40host.com@ftp.ncbi.nih.gov/refseq/release/microbial/"
(defn process-refseq-files [url loc]
  (let [ftp-files (ftp-dir url refseq-files)
        [gbffs fnas] (let [x (group-by #(if (re-find #"\.gbff" %) :gbff :fna)
                                       ftp-files)]
                       [(x :gbff) (x :fna)])
        loc (fs/fullpath loc)
        gbff-dir (str loc "/GBank")
        fna-dir  (str loc "/Fasta")]
    (when (not (fs/directory? gbff-dir))
      (fs/mkdir gbff-dir))
    (when (not (fs/directory? fna-dir))
      (fs/mkdir fna-dir))
    ;; (log> :info "Start downloads")
    (retrieve-refseq-files url gbffs gbff-dir)
    (retrieve-refseq-files url fnas fna-dir)
    ;; (log> :info "Finished downloads")
    :success))


(def refseq-ftp-url
  "ftp://anonymous:user%40host.com@ftp.ncbi.nih.gov/refseq/release/microbial/")

;;; (def *ref-ftp-load*
;;;      (future (process-refseq-files refseq-ftp-url "/data2/BioData")))


(def qin-et-al-2010-url
     "http://www.bork.embl.de/~arumugam/Qin_et_al_2010/")

(defn fetch-qin-fastas [url local-dir]
  (let [toc (fetch-url url)
        files (filter #(re-find #"\.seq\.fa\.gz$" %)
                      (fetch-elements-by-tag toc [:a]))]
    (doseq [f files]
      (http-binary-download
       (str qin-et-al-2010-url f)
       (fs/join local-dir f)))))

;;;(fetch-qin-fastas qin-et-al-2010-url "/data2/Bio/MetaG1/FastaFiles")))




(def rfam-url-template
     ["http://rfam.sanger.ac.uk/family/alignment/download/format?acc="
                      ; RFAM accession
      "&alnType="     ; seed or full
      "&nseLabels="   ; 1 => name/start-end, 0 => species labels
      "&format="      ; stockholm, pfam, fasta, fastau
      "&download=0"])

(defn download-rfam-stos
  "Download the RFAM alignment alignments for the family denoted by
   the RFAM accessions in ACCESSIONS (a seq).  Files are named after
   the accessions and placed in DIR, a filespec for a directory.

   TYPE is either SEED or FULL.  LABELS is 1 for name/start-end or 0
   for species labels.  FORMAT is one of stockholm, pfam, fasta, or
   fastau.  The last, fastau is for ungapped fasta.  Defaults are
   type/seed, labels/1, and format/stockholm"

  [accessions dir
   & {type :type labels :labels format :format
      :or {type "seed" labels 1 format "stockholm"}}]

  (let [origin "#=GF AU Infernal 1.0.2"
        dir (fs/fullpath dir)
        suffixes {"stockholm" ".sto"
                  "pfam" ".pfam"
                  "fasta" ".fna"
                  "fastau" ".fnau"}
        files (map #(fs/join dir (str % "-" type (suffixes format)))
                   accessions)
        urls (map #(str/join
                    "" (interleave rfam-url-template [% type labels format]))
                  accessions)]
    (map #(do (join-sto-fasta-file % % :origin origin) %)
         (map #(do (http-binary-download %1 %2) %2)
              urls
              files))))

;;; (download-rfam-stos
;;;  (keep #(when (> (count %) 1) (first %))
;;;        (csv/parse-csv (slurp "/data2/Bio/ECRibLeaders/SearchRNAs-50-set.csv")))
;;;  "/data2/Bio/ECRibLeaders/RFAM")




(def embl-url-template
     ["http://www.ebi.ac.uk/Tools/dbfetch/emblfetch?id="
      "&Submit=Go"])

(defn embl-acc-ncbi-taxon-id [embl-acc]
  (let [taxon-field
        (re-find #"taxon:[0-9]+"
                 (with-out-str
                   (with-open [is (io/input-stream
                                   (str (embl-url-template 0)
                                        embl-acc ;"AAND01000066"
                                        (embl-url-template 1)))]
                     (io/copy is *out*))))]
    [embl-acc  (if taxon-field
                 (second (str/split #":" taxon-field))
                 ;; Obsolete acc or other issue...
                 "0")]))


(defn get-rfam-family-taxons [rfam-sto-filespec]
  (let [rfam-sto (fs/fullpath rfam-sto-filespec)
        taxon-file (fs/replace-type rfam-sto ".taxons")
        pairs (map embl-acc-ncbi-taxon-id
                   (keep  #(when (not (or (.startsWith % "//")
                                          (.startsWith % "#=GC")))
                             (first (str/split #"/" %)))
                          (first (sto-GC-and-seq-lines rfam-sto-filespec))))]
    (io/with-out-writer taxon-file
      (doseq [[nm tx] pairs]
        (println (str nm "," tx))))))


(defn gen-rfam-family-taxons [rfam-sto-dir]
  (let [base (fs/fullpath rfam-sto-dir)
        stos (sort (fs/directory-files base ".sto"))]
    (doall (map get-rfam-family-taxons stos))))






