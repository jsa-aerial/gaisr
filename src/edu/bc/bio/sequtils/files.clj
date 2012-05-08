;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                      S E Q U T I L S . F I L E S                         ;;
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

(ns edu.bc.bio.sequtils.files

  "Various bio sequence file format readers, writers, verifiers, and
   manipulators."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clj-shell.shell :as sh]
            [edu.bc.fs :as fs])
  (:use clojure.contrib.math
        edu.bc.utils
	edu.bc.utils.probs-stats
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        ))




;;; ----------------------------------------------------------------------
;;;
;;; Convert Sto and Fasta split sequence format files into conjoined
;;; versions.  Many Sto and Fasta files from various sites come in old
;;; fashioned 80 col mode where sequences are split at 80 column mark.
;;; For aligned files this is even worse as you have groups of
;;; sequences split across lines separated by whole pages of other
;;; sequence (parts).  For example, RFAM alignments.  This group puts
;;; all those back together so that each sequence is on a single line.


(defn sto-GC-and-seq-lines [stofilespec]
  (seq/separate
   #(and (> (count %) 1)
         (or (not (.startsWith % "#"))
             (or (.startsWith % "#=GC SS_cons")
                 (.startsWith % "#=GC RF"))))
   (io/read-lines (fs/fullpath stofilespec))))


(defn join-sto-fasta-lines [infilespec origin]
  (let [[seqcons-lines gc-lines] (sto-GC-and-seq-lines infilespec)
        gc-lines (if (not= origin "")
                   (concat (take 1 gc-lines) [origin] (drop 1 gc-lines))
                   gc-lines)
        recombined-lines (sort-by
                          #(-> % second first)
                          (vec (reduce
                                (fn [m l]
                                  (let [[nm sq]
                                        (cond
                                         ;;splits the line apart and
                                         ;;hopefully creates vector
                                         ;;["#GC SS_cons" structure]
                                         (.startsWith l "#=GC SS_cons")
                                         [(str/join " " (butlast (str/split #"\s+" l)))
                                          (last (str/split #"\s+" l))]

                                         (.startsWith l "#")
                                         (str/split #"\s{2,}+" l)

                                         :else
                                         (str/split #"\s+" l))

                                        prev (get m nm [(gen-uid) ""])]
                                    (assoc m  nm [(first prev)
                                                  (str (second prev) sq)])))
                                {} seqcons-lines)))
        {seq-lines false cons-lines true} (group-by
                                           #(or (.startsWith (first %) "//")
                                                (.startsWith (first %) "#"))
                                           recombined-lines)]
    [gc-lines seq-lines cons-lines]))


(defn join-sto-fasta-file
  "Block/join unblocked sequence lines in a sto or fasta file. For sto
   files ORIGIN is a #=GF line indicating tool origin of file.  For
   example, '#=GF AU Infernal 1.0.2'. Defaults to nothing."

  [in-filespec out-filespec
   & {origin :origin :or {origin ""}}]
  (let [[gc-lines seq-lines cons-lines]
        (join-sto-fasta-lines in-filespec origin)]
    (io/with-out-writer (fs/fullpath out-filespec)
      (doseq [gcl gc-lines] (println gcl))
      (doseq [sl seq-lines]
        (let [[nm [_ sq]] sl]
          (cl-format true "~A~40T~A~%" nm sq)))
      (doseq [cl cons-lines]
        (let [[nm [_ sq]] cl]
          (cl-format true "~A~40T~A~%" nm sq))))))


;;; Convert STO format to ALN format (ClustalW format).  This is
;;; needed by some processors which cannot take a Stockholm alignment
;;; format but need an "equivalent" in ClustalW ALigNment format.
;;;
;;; OK, (9-Feb-2012) some tools seem to need things blocked while
;;; others don't work if they are blocked.  Worse, what counts as
;;; valid Clustal/aln format or not is ill defined with multiple
;;; definitions in the community (e.g., many claim 60 character seqs
;;; per line but others say 50; some claim must be blocked, others say
;;; unblocked is valid).  So, we have two variants.  One which blocks
;;; and a main driver which calls blocked version if requested or just
;;; does simple unblocked itself.
;;;
(defn sto->aln-blocked
  "Convert a stockhom format alignment file into its ClustalW
   equivalent BLOCKED ALN format. Blocking is done in 60 character
   chunks.  STOIN is the filespec for the stockholm format file and
   ALNOUT is the filespec for the resulting conversion (it is
   overwritten if it already exists!)"

  [stoin alnout]
  (let [seq-lines (second (join-sto-fasta-lines stoin ""))
        seq-lines (map (fn [[nm [uid sl]]]
                         [nm [uid (partition-stg
                                   60 (str/replace-re #"\." "-" sl))]])
                       seq-lines)]
    (io/with-out-writer alnout
      (println "CLUSTAL W (1.83) multiple sequence alignment\n")
      (loop [x seq-lines]
        (let [[nm [uid sl]] (first x)]
          (when (not-empty sl)
            (do
              (doseq [[nm [uid sl]] x]
                  (cl-format true "~A~40T~A~%" nm (first sl)))
              (println "")
              (recur (map (fn [[nm [uid sl]]]
                            [nm [uid (rest sl)]])
                          x)))))))
    alnout))

(defn sto->aln
  "Convert a stockhom format alignment file into its ClustalW
   equivalent ALN format.  STOIN is the filespec for the stockholm
   format file and ALNOUT is the filespec for the resulting
   conversion (it is overwritten if it already exists!)

   BLOCKED is a boolean indicating whether the output should be
   blocked (60 chars per chunk).  Default is unblocked."

  [stoin alnout & {blocked :blocked :or {blocked false}}]
  (if blocked
    (sto->aln-blocked stoin alnout)
    (let [seq-lines (filter #(not (or (.startsWith % "//") (re-find #"^#" %)))
                            (first (sto-GC-and-seq-lines stoin)))
          seq-lines (map #(str/replace-re #"\." "-" %) seq-lines)]
      (io/with-out-writer (fs/fullpath alnout)
        (println "CLUSTAL W (1.83) multiple sequence alignment\n")
        (doseq [sl seq-lines]
          (println sl)))
      alnout)))


(defn check-sto
  "Checks a sto file to ensure that there are valid characters being
   used in the sequences consensus structure line. Will print out
   errors in the sto file by sequence number.  Input requires a sto
   file"

  [sto & {printem :printem :or {printem true}}]
  (let [valid-symbols #{"A" "C" "G" "U"
                        "-" "." ":" "_" ","
                        "a" "b" "B"
                        "(" ")" "<" ">"}
        [_ seq-lines cons-lines] (join-sto-fasta-lines sto "")
        [_ [_ cl]] (first cons-lines)
        ;;adds namez to the sequences so that they
        ;;can be identified
        sl (map (fn [[nm [_ sq]]]
                  [nm (.toUpperCase sq)])
                seq-lines)

        ;;checks for valid symbols
        check-char (fn [[n s]]
                     [n (every? #(contains? valid-symbols %)
                                (rest (str/split #"" s)))])
        ;;check for common case of repeat named seqs
        check-double-len (fn [[n s]]
                           [n (>= (count s) (* 2 (count cl)))])
        ;;checks to see if sequences have same
        ;;length as consensus structure
        check-len (fn [[n s]]
                    [n (= (count s) (count cl))])

        chks [["sequence contains invalid character in: "
               (map first (remove #(second %) (map #(check-char %) sl)))]
              ["sequence is repeated - two or more times with same name: "
               (map first (filter #(second %) (map #(check-double-len %) sl)))]
              ["sequence contains invalid length compared to cons-line in: "
               (map first (remove #(second %) (map #(check-len %) sl)))]
              ["consensus structure contains invalid character"
               (if-let [x (second (check-char ["#SS_cons" cl]))]
                 ()
                 ["#SS_cons"])]]]

    (when printem
      (doseq [[msg x] chks]
        (when (seq x)
          (prn msg x))))

    (if (every? #(-> % second seq not) chks)
      (do (when printem (prn "sto is good"))
          :good)
      (filter #(seq (second %)) chks))))




;;; ----------------------------------------------------------------------
;;; BLAST oriented file functions.  These are currently based on NCBI
;;; BLAST+.  For BLAST proper support, see sequtils.tools

(def default-binary-db
     (or (getenv "MLAB_DEFAULT_BINARY_DB")
         "/data2/BioData/BlastDBs/other_genomic"))


(defn gen-entry-file [entries file]
  (io/with-out-writer (io/file-str file)
    (doseq [e entries]
      (println e)))
  file)


(defn has-loc? [entries]
  (let [entries (cond
                 (seq? entries) entries
                 (fs/exists? entries) (str/split #"\n" (slurp entries))
                 :else (raise :type :unknown-entries :entries entries))]
    (let [e (first (ensure-vec entries))]
      (re-find #"[0-9]+\-[0-9]+" e))))




(defn nms-sqs->fasta-file [nms-sqs filespec]
  (let [filespec (fs/fullpath filespec)]
    (io/with-out-writer filespec
      (doseq [[nm sq] nms-sqs]
        (println (str ">" nm))
        (println sq)))
    filespec))


(defn entry-range->fasta-file [entry range blastdb & [filespec]]
  (let [fasta-filespec (if filespec filespec (fs/tempfile entry ".fna"))
        blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")]
    (assert-tools-exist [blastdbcmd])
    (let [[range strand] (str/split #"/" range)
          strand (if (= "-1" strand) "minus" "plus") ; => default plus
          cmdargs ["-db" blastdb
                   "-entry" entry "-range" range "-strand" strand
                   "-line_length" "14000000"
                   "-outfmt" "%f" "-out" (fs/fullpath fasta-filespec)]]
      (runx blastdbcmd cmdargs)
      fasta-filespec)))


(defn- entry-file->fasta-file-ranges [efile fasta-filespec blastdb]
  (let [blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")
        tmp-file (fs/tempfile "fasta-out" ".fna")]
    (assert-tools-exist [blastdbcmd])
    (io/with-out-writer fasta-filespec
      (do-text-file [efile]
        (let [[entry range] (str/split #" " $line)
              [range strand] (str/split #"/" range)
              strand (if (= "-1" strand) "minus" "plus") ; => default plus
              cmdargs ["-db" blastdb
                       "-entry" entry "-range" range "-strand" strand
                       "-line_length" "14000000"
                       "-outfmt" "%f" "-out" tmp-file]]
          (catch-all (runx blastdbcmd cmdargs))
          (do-text-file [tmp-file]
                        (println $line)))))
    (fs/rm tmp-file)
    fasta-filespec))


(defn- entry-file->fasta-file-full [efile fasta-filespec blastdb]
  (let [blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")
        cmdargs ["-db" blastdb "-entry_batch" efile
                 "-line_length" "14000000"
                 "-outfmt" "%f" "-out" fasta-filespec]]
    (assert-tools-exist [blastdbcmd])
    (catch-all (runx blastdbcmd cmdargs))
    fasta-filespec))


(defnk entry-file->fasta-file [efile :loc nil :blastdb nil]
  (let [efile (fs/fullpath efile)
        filespec (fs/fullpath (fs/replace-type efile ".fna"))
        blastdb (if blastdb blastdb default-binary-db)]
    (if loc
      (entry-file->fasta-file-ranges efile filespec blastdb)
      (entry-file->fasta-file-full efile filespec blastdb))
    filespec))


(defnk entry-file->blastdb [efile :out nil :blastdb nil :loc nil :ids nil]
  (let [entries-fasta (entry-file->fasta-file efile :loc loc :blastdb blastdb)
        blast-db-file (fs/fullpath
                       (if out out (fs/replace-type efile ".blastdb")))
        blast-path (get-tool-path :ncbi)
        mkdbcmd (str blast-path "makeblastdb")
        cmdargs ["-parse_seqids"
                 "-in" entries-fasta "-dbtype" "nucl"
                 "-out" blast-db-file]
        cmdargs (if ids cmdargs (vec (rest cmdargs)))]
    (assert-tools-exist [mkdbcmd])
    (runx mkdbcmd cmdargs)
    blast-db-file))


(defn get-selection-fna
  "Get a fasta file for SELECTIONS.  If selections is a file, assumes
   it is in fact the fasta file.  If selections is a collection,
   assumes the collection is a set of pairs [nm sq], and converts to
   corresponding fasta file.  Returns full path of result file."
  [selections]
  (if (coll? selections)
    (nms-sqs->fasta-file selections)
    (fs/fullpath selections)))



;;;------------------------------------------------------------------------;;;
;;; Various sequence file readers for various formats.  In particular,
;;; fna, sto, and aln


(defn gaisr-seq-set?
  "Returns true if either

   X is a filespec string with type extension fna, aln, sto, or gma
   X is a java.io.File (presumed to be of one of the above formats)
   X is a collection (presumed to have seqences as elements)
  "
  [x]
  (or (and (string? x)
           (in (fs/ftype x) ["fna" "aln" "sto" "gma"]))
      (isa? (type x) java.io.File)
      (coll? x)))


(defn read-seqs
  "Read the sequences in FILESPEC and return set as a seq (Clojure
  seq!).  Filespec can denote either a fna, aln, sto, or gma file
  format file.
  "
  [filespec]
  (let [type (fs/ftype filespec)
        sqs (drop-until #(re-find #"^(>N|N)" %)
                        (str/split-lines (slurp filespec)))
        sqs (case type
              "aln"
              (map #(str/replace-re #"^N[CZ_0-9]+ +" "" %) sqs)

              "sto"
              (map #(str/replace-re #"^N[CZ_0-9]+ +" "" %)
                   (take-while #(re-find #"^N" %) sqs))

              "gma" (raise :type :NYI :info "GMA format not yet implemented")

              "fna"
              (map second (partition 2 sqs)))]
    sqs))


(defn read-aln-seqs
  "Read the _aligned_ sequences in FILESPEC and return them in a
  Clojure seq.  Filespec can denote either an aln, sto, or gma file
  format file.  If COLS is true, return the _columns_ of the
  alignment (including gap characters).
  "
  [filespec & {cols :cols :or {cols false}}]
  {:pre [(in (fs/ftype filespec) ["aln" "sto" "gma"])]}
  (let [sqs (read-seqs filespec)]
    (if cols (transpose sqs) sqs)))

(defn read-dirs-aln-seqs
  "Apply read-aln-seqs across FTYPES in DIR.  FTYPES is a vector of
   one or more file formats (as type extensions) taken from #{\"fna\",
   \"aln\", \"sto\", \"gma\"}, producing a seq of sequence sets, each
   set being a seq of the sequences in a matching file.  NOTE: each
   such set can be viewed as a 'matrix' of the elements (bases/gaps)
   of the sequences.

   If COLS is true, return the transpose of each obtained sequence
   'matrix', i.e., return the seq of column seqs.  NOTE: as the
   formats are alignment formats, all sequences in a file are of the
   same length (including gaps).

   If DIRDIR is true, dir is taken as a directory of directories, each
   of which will have read-aln-seqs applied to it per the above
   description and the return will be a seq of all such applications.

   So, if dirdir is false, the result will have the form:

   (seqs-from-file1 ... seqs-from-filen), where filei is in dir

   If dirdir is true, the result will be nested one level more:

   ((seqs-from-dir1-file1 ... seqs-from-dir1-filek)
    ...
    (seqs-from-dirn-file1 ... seqs-from-dirn-filel))
  "
  [dir & {dirdir :dirdir cols :cols ftypes :ftypes
          :or {dirdir false cols false ftypes ["aln"]}}]
  (let [one-dir (fn[dir]
                  (filter #(> (count %) 1)
                          (fs/dodir
                           dir (fn[d]
                                 (apply set/union
                                        (map #(fs/directory-files d %)
                                             ftypes)))
                           read-aln-seqs)))
        sqs (if dirdir
              (map one-dir (sort (fs/directory-files dir "")))
              (one-dir dir))]
    (if cols
      (if dirdir
        (map #(map transpose %) sqs)
        (map transpose sqs))
      sqs)))
