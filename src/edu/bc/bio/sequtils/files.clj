;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                      S E Q U T I L S . F I L E S                         ;;
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
;; Author: Jon Anthony, Shermin Pei                                         ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns edu.bc.bio.sequtils.files

  "Various bio sequence file format readers, writers, verifiers, and
   manipulators."

  (:require [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.contrib.io :as io]
            [clojure-csv.core :as csv]
            [edu.bc.fs :as fs])
  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        [edu.bc.bio.seq-utils :only [reverse-compliment]]
        [clojure.pprint
         :only [cl-format]]
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
   (filter #(not= (str/replace-re #"\s+" "" %) "")
           (io/read-lines (fs/fullpath stofilespec)))))


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
                                         [(str/join " " (butlast (str/split
                                                                  #"\s+" l)))
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
  "Block/join unblocked sequence lines in a sto file.  ORIGIN is a
   #=GF line indicating tool origin of file.  For example, '#=GF AU
   Infernal 1.0.2'. Defaults to nothing."

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


(def default-genome-fasta-dir
     (or (getenv "GAISR_DEFAULT_GENOME_DIR")
         "/data2/BioData/GenomeSeqs/RefSeq58"))

(defn split-join-fasta-file
  [in-file out-dir
   & {:keys [base pat] :or {base default-genome-fasta-dir pat #"^gi"}}]
  (doseq [[gi sq] (->> in-file io/read-lines
                       (partition-by #(re-find pat %))
                       (partition-all 2)
                       (map (fn[[[nm] sbits]] [nm (apply str sbits)])))]
    (let [nm (->> gi (str/split #"\|") (drop 3) first
                  (str/split #"\.") first)]
      (when (re-find #"^NC_" nm)
        (io/with-out-writer (fs/join base (str nm ".fna"))
          (println gi)
          (println sq))))))

(defn split-join-ncbi-fasta-file
  "Split a fasta file IN-FILE into the individual sequences and
   unblock the sequence if blocked.  The resulting individual [nm sq]
   pairs are written to files named for the NC name in the gi line of
   in-file and in the DEFAULT-GENOME-FASTA-DIR location.

   The main use of this function is to take a refseq fasta
   db (composed of many multi seq fasta files) and split the db into a
   normed set of named sequence files for quick access to sequence per
   name in various other processing (see gen-name-seq for example).

   Canonical use case example:

   (fs/dodir \"/data2/BioData/Fasta\" ; RefSeqxx fasta files
             #(fs/directory-files % \"fna\")
             #(split-join-ncbi-fasta-file %))
  "
  [in-file]
  (let [base default-genome-fasta-dir]
    (doseq [[gi sq] (->> in-file io/read-lines
                         (partition-by #(re-find #"^>gi" %))
                         (partition-all 2)
                         (map (fn[[[nm] sbits]] [nm (apply str sbits)])))]
      (let [nm (->> gi (str/split #"\|") (drop 3) first
                    (str/split #"\.") first)]
        (when (re-find #"^NC_" nm)
          (io/with-out-writer (fs/join base (str nm ".fna"))
            (println gi)
            (println sq)))))))


(defn chunk-genome-fnas
  "Take the set of fnas in directory GENOME-FNA-DIR (presumably
   created by split-join-ncbi-fasta-file or similar) and aggregate
   them into a new smaller set of files, where each new file contains
   the contents of CHUNK-SIZE input files (with the possible exception
   of the last file having a smaller number).  This is useful for
   creating custom data sets for search.
  "
  [genome-fna-dir &
  {:keys [chunk-size] :or {chunk-size 100}}]
  (let [dir (fs/join genome-fna-dir "Chunked")
        all (->> (fs/directory-files genome-fna-dir ".fna") sort
                 (partition-all chunk-size))]
    (when (not (fs/exists? dir)) (fs/mkdir dir))
    (doseq [grp all]
      (let [n1 (-> grp first fs/basename (fs/replace-type ""))
            n2 (-> grp last fs/basename (fs/replace-type ""))
            file (fs/join dir (str n1 "-" n2 ".fna"))]
        (io/with-out-writer file
          (doseq [f grp
                  l (io/read-lines f)]
            (println l)))))))


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




;;; Cool stuff from Shermin, but requires much more refactoring of
;;; various other things of his to make it all work. See his fold-ops.
(defn print-sto
  "takes sequence lines and a structure line and writes it into a sto
  format file. the seq-lines needs to be a collection of [name
  sequence] pairs. structure is a string. Simply prints out to the
  repl."

  [seq-lines structure]
  (println "# STOCKHOLM 1.0\n")
  (doseq [sq seq-lines]
    (let [[nm sq] (if (vector? sq)
                    sq
                    (str/split #"\s+" sq))]
      (cl-format true "~A~40T~A~%" nm (str/replace-re #"\-" "." sq))))
  (cl-format true "~A~40T~A~%" "#=GC SS_cons" structure)
  (println "//"))

#_(defn aln->sto
  "takes an alignment in Clustal W format and produces a sto file by
   using RNAalifold to determine the structure and then making it into
   a sto file adding header and a consensus line"

  [in-aln out-sto & {:keys [fold-alg st]
                     :or {fold-alg "RNAalifold"}}]
  (cond
   (identity st) ;structure provided
   (let [sq (read-seqs in-aln :type "aln")]
     (io/with-out-writer out-sto
       (print-sto sq st))
     out-sto)

   (= fold-alg "RNAalifold") ;structure from RNAalifold
   (let [st (fold-aln in-aln)
         sq (read-seqs in-aln :type "aln")]
     (io/with-out-writer out-sto
       (print-sto sq st))
     out-sto) ;return out sto filename

   ;;else use cmfinder
   #_(shell/sh "perl" "/home/kitia/bin/gaisr/src/mod_cmfinder.pl" in-aln out-sto)))



;;; Forward declarations...
(declare
 read-seqs
 entry-parts)

(defn sto->fna
  "Convert a sto file into a fasta file.  Split seq lines into names
   and seq data and interleave these.  Seq data has all gap characters
   removed."
  [stoin fnaout]
  (io/with-out-writer fnaout
    (doseq [[nm seqdata] (read-seqs stoin :info :both)]
      (println (str ">" nm))
      (println (str/replace-re #"[-.]+" "" seqdata))))
  fnaout)


(defn check-sto
  "Checks a sto file to ensure that there are valid characters being
   used in the sequences consensus structure line. Will print out
   errors in the sto file by sequence number.  Input requires a sto
   file"

  [sto & {printem :printem :or {printem true}}]
  (let [valid-symbols #{"A" "C" "G" "U"
                        "-" "." ":" "_" ","
                        "a" "b" "B" "n" "N"
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

        ;; Check for valid entry formats.  Defined as parseable by
        ;; entry-parts
        check-entry (fn [[n _]]
                      [n (not (vector? (-> n entry-parts catch-all)))])

        chks [["sequence contains invalid character in: "
               (map first (remove #(second %) (map #(check-char %) sl)))]
              ["sequence is repeated - two or more times with same name: "
               (map first (filter #(second %) (map #(check-double-len %) sl)))]
              ["sequence contains invalid length compared to cons-line in: "
               (map first (remove #(second %) (map #(check-len %) sl)))]
              ["sequence line has invalid entry - not parseable: "
               (map first (remove #(not (second %)) (map #(check-entry %) sl)))]
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
     (or (getenv "GAISR_DEFAULT_BINARY_DB")
         "/data2/BioData/BlastDBs/RefSeq58/refseq58_microbial_genomic"))


(defn make-entry
  "'Inverse' of entry-parts.  EVEC is a vector of shape [nm [s e]
   strd], where nm is the entry name, S and E are the start and end
   coordinates in the genome, and strd is the strand marker, 1 or
   -1. Returns the full entry as: nm/s-e/strd
  "
  ([evec]
     (let [[nm [s e] st] evec]
       (str nm "/" s "-" e "/" st)))
  ([nm s e st]
     (make-entry [nm [s e] st])))


(defn entry-parts
  "ENTRY is a string \"name/range/strand\", where name is a genome
   name, range is of the form start-end and strand is 1 or -1.  At
   least name must be supplied.  DELTA is an integer, which will be
   subtracted from start and added to end.

   Returns a triple [name [start end] strand]
  "
  [entry & {:keys [ldelta rdelta] :or {ldelta 0 rdelta 0}}]
  (let [[name range] (str/split #"( |/|:)+" 2 entry)
        name (re-find #"[A-Za-z0-9]+_[A-Za-z0-9_]+" name)
        [range strand] (if range (str/split #"/" range) [nil nil])
        [s e st] (if (not range)
                   [1 Long/MAX_VALUE "1"]
                   (let [range (if (= \c (.charAt range 0))
                                 (subs range 1)
                                 range)
                         [s e] (map #(Integer. %) (str/split #"-" range))
                         strand (if strand strand (if (> s e) "-1" "1"))
                         [s e] (if (< s e) [s e] [e s])
                         [s e] [(- s ldelta) (+ e rdelta)]
                         [s e] [(if (<= s 0) 1 s) e]]
                     [s e strand]))]
    [name [s e] st]))


(defn gen-entry-file [entries file]
  (io/with-out-writer (io/file-str file)
    (doseq [e entries]
      (println e)))
  file)

(defn gen-entry-nv-file [entries file]
  (io/with-out-writer (io/file-str file)
    (doseq [e entries]
      (if (coll? e)
        (println (->> e (map str) (str/join ", ")))
        (println e))))
  file)


(defn has-loc? [entries]
  (let [entries (cond
                 (seq? entries) entries
                 (fs/exists? entries) (str/split #"\n" (slurp entries))
                 :else (raise :type :unknown-entries :entries entries))]
    (let [e (first (ensure-vec entries))]
      (re-find #"[0-9]+\-[0-9]+" e))))




(defn gen-name-seq
  "Generate a pair [entry genome-seq], from ENTRY as possibly modified
   by [L|R]DELTA and RNA.  ENTRY is a string \"name/range/strand\",
   where

   name is the genome NC name (and we only currently support NCs),

   range is of the form start-end, where start and end are integers (1
   based...) for the start and end coordinates of name's sequence to
   return.  NOTE: start < end, as reverse compliment information comes
   from strand.

   strand is either -1 for reverse compliment or 1 for standard 5'->3'

   LDELTA is an integer, which will be _subtracted_ from start.  So,
   ldelta < 0 _removes_ |ldelta| bases from 5', ldelta > 0 'tacks on'
   ldelta extra bases to the 5' end.  Defaults to 0 (no change).

   RDELTA is an integer, which will be _added_ to the end.  So, rdelta
   < 0 _removes_ |rdelta| bases from 3', rdelta > 0 'tacks on' rdelta
   extra bases to the 3' end.  Defaults to 0 (no change).

   If RNA is true, change Ts to Us, otherwise return unmodified
   sequence.

   BASEDIR is the location of the NC fasta files.  Generally, this
   should always be the default location.
  "
  [entry & {:keys [basedir ldelta rdelta rna]
            :or {basedir default-genome-fasta-dir ldelta 0 rdelta 0 rna true}}]
  (let [[name [s e] strand] (entry-parts entry :ldelta ldelta :rdelta rdelta)
        _ (when (<= s 0)
            (raise :type :bad-entry
                   :input [entry ldelta rdelta]
                   :result [name s e strand]))
        entry (str name "/" s "-" e "/" strand)
        fname (fs/join basedir (str name ".fna"))
        sq (->> (io/read-lines fname) second
                (str/drop (dec s)) (str/take (- (inc e) s)))
        sq (if (= strand "-1") (reverse-compliment sq) sq)
        sq (if rna (str/replace-re #"T" "U" sq) sq)]
    [entry sq]))


;;; (fs/dodir
;;;  "/home/kaila/Bio/Work/S15work/062612"
;;;  #(fs/directory-files % ".ent")
;;;  gen-name-seq-pairs)
;;;
(defn gen-name-seq-pairs
  "Generate a sequence of pairs [name genome-seq] from the given
   collection denoted by entries, which is either a collection or a
   string denoting an entry filespec (either hand made or via
   gen-entry-file or similar).

   See gen-name-seq for details.  Basically this is (map gen-name-seq
   entries) with some options.  In particular, entries may lack a
   strand component (again see gen-name-seq for format details), and
   STRAND here would be used to indicate all entries are to have this
   strand.  If entries have a strand component, strand should be 0,
   otherwise you will force all entries to the strand given.  Strand
   is either -1 (reverse compliment) or 1 for standard 5'->3', or 0
   meaning 'use strand component in entries'

   BASEDIR is the location of the NC fasta files.  Generally, this
   should always be the default location.
  "
  [entries & {:keys [basedir strand ldelta rdelta rna]
              :or {basedir default-genome-fasta-dir
                   strand 0 ldelta 0 rdelta 0 rna true}}]
  (let [entries (if (string? entries) (io/read-lines entries) entries)
        entries (if (= 0 strand) entries (map #(str % "/" strand) entries))]
    (map #(gen-name-seq
           % :basedir basedir :ldelta ldelta :rdelta rdelta :rna rna)
         entries)))




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


(defn- ncbi-entry-file->fasta-file-ranges
  "Use NCBI blastdbcmd to get sequences from coordinates - OBSOLETE!"
  [efile fasta-filespec blastdb]
  (let [blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")
        tmp-file (fs/tempfile "fasta-out" ".fna")]
    (assert-tools-exist [blastdbcmd])
    (io/with-out-writer fasta-filespec
      (do-text-file [efile]
        (let [[entry range] (str/split #"( |/)+" 2 $line)
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

(defn- ncbi-entry-file->fasta-file-full
  "Use NCBI blastdbcmd to get full sequences - OBSOLETE!"
  [efile fasta-filespec blastdb]
  (let [blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")
        cmdargs ["-db" blastdb "-entry_batch" efile
                 "-line_length" "14000000"
                 "-outfmt" "%f" "-out" fasta-filespec]]
    (assert-tools-exist [blastdbcmd])
    (catch-all (runx blastdbcmd cmdargs))
    fasta-filespec))

(defn ncbi-entry-file->fasta-file
  "Use NCBI blastdbcmd to get fasta file for entries in entry file - OBSOLETE"
  [efile & {:keys [loc blastdb] :or {loc nil blastdb default-binary-db}}]
  (let [efile (fs/fullpath efile)
        filespec (fs/fullpath (fs/replace-type efile ".fna"))]
    (if loc
      (ncbi-entry-file->fasta-file-ranges efile filespec blastdb)
      (ncbi-entry-file->fasta-file-full efile filespec blastdb))
    filespec))


(defn entry-file->fasta-file
  [efile & {:keys [names-only]}]
  (let [efile (fs/fullpath efile)
        fasta-filespec (fs/fullpath (fs/replace-type efile ".fna"))
        entseqs (if names-only
                  (read-seqs efile :info :both)
                  (gen-name-seq-pairs (io/read-lines efile)))]
    (io/with-out-writer fasta-filespec
      (doseq [[entry sq] entseqs]
        (println (str ">" entry))
        (println sq)))
    fasta-filespec))


(defn entry-file->blastdb
  [efile & {:keys [out ids] :or {out nil ids nil}}]
  (let [entries-fasta (entry-file->fasta-file efile)
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


(defn gen-nc-genome-fnas
  "One shot genome sequence fnas generator.  Typically used once per
   data update.  Needs to be generalized to be able to use Genbank fna
   archives.  Currently assumes one fna AND ONLY NC_* genomes.
  "
  [full-nc-genomes-filespec]
  (let [dir (fs/dirname full-nc-genomes-filespec)
        fna-pairs (partition-all 2 (io/read-lines full-nc-genomes-filespec))]
    (doseq [[gi genome] fna-pairs]
      (let [fname (fs/join dir (str (re-find #"NC_[0-9]+" gi) ".fna"))]
        (io/with-out-writer fname
          (println gi)
          (println genome))))))




;;; ----------------------------------------------------------------------- ;;;

;;;  This crap needs to be refactored and most of it probably
;;;  eliminated.  Most of it is old stuff that is largely superceded
;;;  but still required by things in post-db-csv, but we at least fix
;;;  up the names of the csv processors to more accurately reflect
;;;  what they are!


(defn get-legacy-csv-info [rows]
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

(defn get-gaisr-csv-info [rows]
  (loop [entries []
         rows (drop 1 rows)]
    (let [entry (first rows)]
      (if (< (count entry) 12) ; minimum used fields
        (sort #(string-less? (%1 0) (%2 0)) entries)
        (recur (conj entries
                     [(entry 0) (entry 3) (entry 4) (entry 9)
                      (entry 8)
                      :new "N/A"])
               (drop 1 rows))))))


(defn canonical-csv-entry-info [entries]
  (map #(let [[nm [s e] sd] (entry-parts %)
              [s e] (if (= sd "1") [s e] [e s])]
          [nm s e 0.0 0.0 :new sd])
       entries))

(defn get-sto-as-csv-info [stofile]
  (let [entries (read-seqs stofile :info :name)]
    (canonical-csv-entry-info entries)))

(defn get-ent-as-csv-info [ent-file]
  (let [entries (io/read-lines ent-file)]
    (if (> (->> entries first csv/parse-csv first count) 1)
      (let [einfo (canonical-csv-entry-info
                   (map #(->> % (str/split #",") first) entries))
            entropy-scores (->> ent-file slurp csv/parse-csv
                                butlast (map second))]
        (map (fn[[nm s e _ _ x y] score] [nm s e score 0.0 x y])
             einfo entropy-scores))
      (canonical-csv-entry-info entries))))


(defn get-csv-entry-info [csv-hit-file]
  (let [file (fs/fullpath csv-hit-file)
        ftype (fs/ftype file)]
    (cond
     (= ftype "sto") (get-sto-as-csv-info file)
     (= ftype "ent") (get-ent-as-csv-info file)
     :else
     (let [rows (csv/parse-csv (slurp file))
           head (first rows)]
       (if (= (first head) "gaisr name")
         (get-gaisr-csv-info rows)
         (get-legacy-csv-info rows))))))



;;; This one was from dists and punted to the old get-enties in
;;; post-db-csv for csv entries
;;;
(defn get-entries
  [filespec & [seqs]]
  (let [fspec (fs/fullpath filespec)
        ftype (fs/ftype fspec)]
    (if (not= ftype "csv")
      (read-seqs filespec :info (if seqs :data :name))
      (->> (get-csv-entry-info fspec)
           (keep (fn[[nm s e sd]]
                   (when (fs/exists? (fs/join default-genome-fasta-dir
                                              (str nm ".fna")))
                     (str nm "/"
                          (if (> (Integer. s) (Integer. e))
                            (str e "-" s "/-1")
                            (str s "-" e "/1"))))))
           (#(if seqs (map second (gen-name-seq-pairs %)) %))))))




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


(defn seqline-info-mapper
  "Helper function for READ-SEQS.  Returns the function to map over
   seq lines to obtain the requested info.  TYPE is supported seq file
   type (aln, sto, fna, fa, gma).  INFO is either :name for the
   sequence identifier, :data for the sequence data, or :both for name
   and data.

   Impl Note: while this almost begs for multimethods, that would
   actually increase the complexity as it would mean 14 methods to
   cover the cases...
  "
  [type info]
  (if (= info :both)
    #(do [((seqline-info-mapper type :name) %)
          ((seqline-info-mapper type :data) %)])
    (case type
          "aln"
          (if (= info :data)
            #(str/replace-re #"^N[CZ_0-9]+\s+" "" %)
            #(re-find  #"^N[CZ_0-9]+" %))

          "sto"
          (if (= info :data)
            #(str/replace-re #"^(N[CZ_0-9]+|[A-Za-z0-9._/-]+)[,\s]+" "" %)
            #(second (re-find  #"^(N[CZ_0-9]+|[A-Za-z0-9._/-]+)[,\s]+" %)))

          "ent"
          ;; This is a bit annoying as we need to account for ent
          ;; files with csv annotation beyond entries
          (if (= info :data)
            #(-> % csv/parse-csv ffirst gen-name-seq second)
            #(->> (str/split #"([\s,]|/)" %) (take 3) (str/join "/")))

          "gma" (raise :type :NYI :info "GMA format not yet implemented")

          ("fna" "fa" "hitfna")
          (if (= info :data)
            second
            #(re-find #"[A-Za-z0-9._/-]+" (first %))))))

(defn read-seqs
  "Read the sequences in FILESPEC and return set as a lazy (Clojure!)
   seq.  Filespec can denote either a fna, fa, hitfna, aln, sto, or
   gma file format file.
  "
  [filespec & {:keys [info type] :or {info :data type (fs/ftype filespec)}}]
  (when (not (fs/empty? filespec))
    (let [f   (seqline-info-mapper type info)
          sqs (filter #(let [l (str/replace-re #"^\s+" "" %)]
                         (and (not= l "") (not (.startsWith l "#"))))
                      (io/read-lines filespec))
          sqs (if (re-find #"^CLUSTAL" (first sqs)) (rest sqs) sqs)
          sqs (drop-until #(re-find #"^(>[A-Za-z]|[A-Za-z])" %) sqs)
          sqs (if (in type ["fna" "fa" "hitfna"]) (partition 2 sqs) sqs)
          sqs (if (= type "sto") (take-while #(re-find #"^[A-Z]" %) sqs) sqs)]
      (map f sqs))))


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
          :or {dirdir false cols false ftypes ["sto"]}}]
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


(defn map-aln-seqs
  ""
  ([f cols filespec]
     (list (f (read-aln-seqs filespec :cols cols))))
  ([f par cols filespec & filespecs]
     (pxmap f par (map #(read-aln-seqs % :cols cols)
                       (cons filespec filespecs)))))

(defn reduce-aln-seqs
  ""
  ([f fr cols filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr (apply map-aln-seqs f par cols filespecs))))
  ([f fr v cols filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr v (apply map-aln-seqs f par cols filespecs)))))


(defn map-seqs
  ([f filespec]
     (list (f (read-seqs filespec))))
  ([f par filespec & filespecs]
     ;;(prn :par par :filespec filespec :filespecs filespecs)
     (pxmap f par (map #(read-seqs %)
                       (cons filespec filespecs)))))

(defn reduce-seqs
  ""
  ([f fr filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr (apply map-seqs f par filespecs))))
  ([f fr v filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr v (apply map-seqs f par filespecs)))))


