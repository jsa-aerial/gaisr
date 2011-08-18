(ns edu.bc.bio.readers-genbank
  "Readers of GBank format files and loaders from these to DBs"
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.xml :as xml]
            [clojure.java.shell :as sh]
            [edu.bc.fs :as fs])
  (:use edu.bc.utils)
  (:import [java.util StringTokenizer]
           [java.io StreamTokenizer]))


(def gbf (io/file-str (fs/fullpath "/data2/BioData/GBank/jsa-tiny.gbff")))

(defn get-spaces [scanner]
  (loop [cnt 1]
    (if (= 32 (. scanner nextToken))
      (recur (inc cnt))
      (do (. scanner pushBack) (str/repeat cnt " ")))))


(defn make-scanner [filespec-or-reader]
  (let [fr filespec-or-reader
        r (if (or (isa? (class fr) java.io.InputStream)
                  (isa? (class fr) java.io.Reader))
            fr
            (io/input-stream (io/file-str (fs/fullpath fr))))
        scanner (StreamTokenizer. r)]
    (doto scanner
      .resetSyntax
      (.eolIsSignificant true)
      (.wordChars (int \_) (int \_))
      (.wordChars (int \A) (int \Z))
      (.wordChars (int \a) (int \z))
      (.wordChars (int \0) (int \9))
      (.wordChars (int \-) (int \-))
      (.wordChars (int \u00A0) (int \u00FF))
      (.quoteChar (int \"))
      (.quoteChar (int \'))
      (.ordinaryChar (int \ ))
      (.ordinaryChar (int \.)))
    scanner))



(defn scan [scanner]
  (let [tk (. scanner nextToken)]
    (cond
     (= tk StreamTokenizer/TT_WORD) [:word (. scanner sval)]
     (= tk StreamTokenizer/TT_NUMBER) [:num (. scanner nval)]
     (= tk StreamTokenizer/TT_EOF) [:eof nil]
     (= tk StreamTokenizer/TT_EOL) [:eol nil]
     :otherwise
     (let [chtk (char tk)]
       (case chtk
             \" [:dq (. scanner sval)]
             \' [:sq (. scanner sval)]
             \space [:spaces (get-spaces scanner)]
             \. [:dot nil]
             \: [:colon nil]
             [(keyword (str chtk)) :other])))))


(defn count-words-lines [filespec]
  (with-open [fd (io/input-stream (io/file-str (fs/fullpath filespec)))]
    (let [scanner (make-scanner fd)]
      (loop [wordcnt 0
             numbercnt 0
             linecnt 0
             [tk v] (scan scanner)]
        (if (not= tk :eof)
          (recur (if (= tk :word) (inc wordcnt) wordcnt)
                 (if (= tk :num) (inc numbercnt) numbercnt)
                 (if (= tk :eol) (inc linecnt) linecnt)
                 (scan scanner))
          {:words wordcnt :numbers numbercnt :lines linecnt})))))



(with-in-str "LOCUS       NC_010628            8234322 bp    DNA     circular BCT 14-MAY-2010.\n       .\n        \"this and that at the 5' end of seq\""
  (let [scanner (make-scanner *in*)]
    (loop [[tk v] (scan scanner)]
      (prn (if (= v :other)
             (str "Ordinary char:" tk)
             (if (= tk :spaces)
               [(str tk ":" (count v)) v]
               (if (not= v nil)
                 (str tk ": " v)
                 (str tk)))))
      (if (not= tk :eof)
        (recur (scan scanner))))))


(with-open [r (io/input-stream gbf)]
  (let [scanner (make-scanner r)]
    (loop [wordcnt 0
           numbercnt 0
           linecnt 0
           [tk v] (scan scanner)]
      (if (= tk :spaces)
        (println (str "Spaces:" (count v) " '" v "'"))
        (if (not= v nil)
          (println tk ":" v)
          (println tk)))
      (if (not= tk :eof)
        (recur (if (= tk :word) (inc wordcnt) wordcnt)
               (if (= tk :num) (inc numbercnt) numbercnt)
               (if (= tk :eol) (inc linecnt) linecnt)
               (scan scanner))
        (str "Words: " wordcnt
             " and Numbers: " numbercnt
             " and Lines: " linecnt)))))

(def cnt 0)
(defn read-genbank [filespec]
  (binding [cnt 0]
    (do-text-file
     [filespec]
     (set! cnt 0)
     (let [bits (str/split #" +" $line)]
       (prn bits)))))

(read-genbank (fs/fullpath "/data2/BioData/GBank/jsa.genomic.gbff"))

