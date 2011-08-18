;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                          D B - D O W N L O A D S                         ;;
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
;; Author                                                         ;;
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
            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.net-utils
        [edu.bc.log4clj :only [create-loggers log>]]

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
