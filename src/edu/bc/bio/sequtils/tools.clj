;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                     S E Q U T I L S . T O O L S                          ;;
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

(ns edu.bc.bio.sequtils.tools

  "Intended to include/contain various utilities for operating with
   and on sequence data.  Originally this was all about working
   with (NCBI) Blast+, CMFINDER (and associated utils) and Infernal
   tools, but now may branch out to many others."

  (:require [clojure.set :as set]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure-csv.core :as csv]
            [edu.bc.fs :as fs])
  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        [edu.bc.bio.seq-utils :only [norm-elements degap-seqs]]
        edu.bc.bio.sequtils.files
        [clojure.pprint
         :only [cl-format]]
        ))


;;; ----------------------------------------------------------------------
;;; NCBI+ BLAST tool access.


(defnk blast
  "Implementation of blast requests.  See blastn, tblastn, etc. for details
   on main arguments.  PGM is a keyword denoting which blast program."
  [pgm in
   :out nil
   :blastdb nil
   :strand :plus
   :word-size 8
   :evalue 10
   :fmt "10 qseqid qstart qend evalue sseqid sstart send"
   :misc nil
   :fn nil
   :fnargs nil]
  (let [seqin       (fs/fullpath in)
        blast-out   (or (and out (fs/fullpath out))
                        (str seqin ".blast"))
        seq-blastdb (or (and blastdb (fs/fullpath blastdb))
                        default-binary-db)

        blast-path  (get-tool-path :ncbi)
        [blastn
         tblastn]   (map #(str blast-path %) '("blastn" "tblastn"))

        io-args     ["-query" seqin "-db" seq-blastdb "-out" blast-out]
        ctrl-args   ["-strand" (name strand)
                     "-word_size" (str word-size)
                     "-evalue" (str evalue)]
        fmt-args    ["-outfmt" fmt]
        misc-args   (if misc (str/split #" " misc) [])]

    (assert-tools-exist [blastn tblastn])

    (case pgm
     :blastn
     (runx blastn (concat io-args ctrl-args fmt-args misc-args))

     :tblastn
     (let [ctrl-args (vec (drop 2 ctrl-args))]
       (runx tblastn (concat io-args ctrl-args fmt-args misc-args))))

    (if-not fn
      blast-out
      (do
        (assert (fn? fn))
        (let [fnargs (if (sequential? fnargs) fnargs (list fnargs))]
          (apply fn blast-out (seq fnargs)))))))


(defn blastn
  "NCBI blastn (from the blast+ toolset).

   IN is the input filespec in fasta format.

   KEYWORD ARGS:
   :out is the output filespec for direct bast output, default `in`.blast
   :blastdb the filespec for the preformatted blast database to use
   :strand (string) plus or minus, defaults to plus
   :word-size blast bp slice size (chunk size for alignment), default 8
   :evalue (real) expectation (E) value threshold for saving hits, default 10
   :fmt NCBI blastn out fmt string, defaults to:
        10 qseqid qstart qend evalue sseqid sstart send
   :fn post processing fn of at least 1 arg for blast output, eg parse-blast
   :fnargs any other args for fn, e.g., output filespec"
  [in & kw-args]
  (apply blast :blastn in kw-args))


(defn tblastn
  "NCBI tblastn (from the blast+ toolset).

   IN is the input filespec in fasta format.

   KEYWORD ARGS:
   :out is the output filespec for direct blast output, default `in`.blast
   :blastdb the filespec for the preformatted blast database to use
   :word-size blast bp slice size (chunk size for alignment), default 8
   :evalue (real) expectation (E) value threshold for saving hits, default 10
   :fmt NCBI blastn out fmt string, defaults to:
        10 qseqid qstart qend evalue sseqid sstart send
   :fn post processing fn of at least 1 arg for blast output, eg parse-blast
   :fnargs any other args for fn, e.g., output filespec"
  [in & kw-args]
  (apply blast :tblastn in kw-args))


(defn blastpgm
  "Determine blast program to use based on sequence alphabet of FASTA-FILE.
   If amino acid alphabet, use tblastn, if nucleotide use blastn.  Returns
   the function for one of these."
  [fasta-file]
  (with-open [r (io/reader fasta-file)]
    (let [l (doall (second (line-seq r)))]
      (if (some #(not (in (first (clojure.string/lower-case %)) "atgcu")) l)
        tblastn
        blastn))))




;;; ----------------------------------------------------------------------
;;; Miscellaneous - I suppose...


;;; From time to time we end up with (or acquire from others) stos
;;; whose seqs have no coordinate info or broken coordinate info.
;;; This takes such stos and regens them with coordinates matching the
;;; seqs.
;;;

(defn genome-seq-dup?
  [coord1 coord2 & {:keys [cutoff] :or {cutoff 0.8}}]

  (let [coord-pat #"([0-9]+)-([0-9]+)"
        [s1 e1] (map #(Integer. %) (rest (re-find coord-pat coord1)))
        sd1 (last (str/split #"/" coord1))
        [s2 e2] (map #(Integer. %) (rest (re-find coord-pat coord2)))
        sd2 (last (str/split #"/" coord2))]
    (cond
     (not= sd1 sd2) ; strands different
     false

     (or (< s1 s2 e2 e1) (< s2 s1 e1 e2))  ; total overlap
     true

     (or (< e1 s2) (< e2 s1)) ; 0 overlap
     false

     :else
     (let [l1 (inc (- e1 s1))
           l2 (inc (- e2 s2))
           lavg (/ (+ l1 l2) 2.0)
           overlap (inc (if (<= s1 s2) (- e1 s2) (- e2 s1)))]
       [l1 l2 lavg overlap (/ overlap lavg)
        (not (neg? overlap))
        (>= (/ overlap lavg) cutoff)]))))

(defn get-blast-id-info-new
  [blast-output-file]
  (reduce
   (fn[m x]
     (let [nm (re-find #"^[A-Za-z_0-9]+" (first x))
           ev (Double. (nth x 3))]
       (if (and (= (second x) "1")
                (str/substring? nm (nth x 4))
                ;; somewhat arbitrary, but real examples should
                ;; have very low evalues
                (< ev 1e-05))
         (let [s (Integer. (nth x 5))
               e (Integer. (nth x 6))
               [s e sd] (if (> s e) [e s -1] [s e 1])
               entry (str nm "/" s "-" e "/" sd)
               sq (if (not (re-find #"^NC" entry))
                    "Non NC_* unavailable"
                    (second (gen-name-seq entry)))]
           (assoc m nm (conj (get m nm [])
                             [(str "/" s "-" e "/" sd) sq])))
         m)))
   {} (butlast (csv/parse-csv (slurp blast-output-file)))))

(defn filter-base-hits-to-most-likely
  [blast-output-file]
  (into {} ; Place in set where keys are names and val is vec of coord/sq pairs
        (for [[nm info] (get-blast-id-info-new blast-output-file)]
          [nm (vec (let [combos (combins 2 info)]
                     (if combos
                       (reduce (fn [m [[c1 sq1] [c2 sq2]]]
                                 (if (genome-seq-dup? c1 c2)
                                   (let [e1 (get m c1)
                                         e2 (get m c2)]
                                     (if (and e1 e2) (dissoc m c2) m))
                                   ;;else, yes we may overwrite same
                                   ;;entry with same value a few times
                                   (assoc (assoc m c1 sq1) c2 sq2)))
                               {} combos)
                       info)))])))

(defn entry-aln-sqs
  [id-info idseq]
  (let [[good bad]
        (reduce (fn[[good bad] [nm coord aln-sq]]
                  (if-let [info (get id-info nm)]
                    (let [base-sq (norm-elements (degap-seqs aln-sq))
                          match-coord (some #(when (= (second %) base-sq)
                                               (first %))
                                            info)]
                      (if match-coord
                        [(assoc good (str nm match-coord)
                                [match-coord coord aln-sq])
                         bad]
                        [good (assoc bad (str nm coord) aln-sq)]))
                     [good (assoc bad (str nm coord) aln-sq)]))
                [{} {}] idseq)]
    [(sort-by key good) (sort-by key bad)]))


(defn correct-sto-coordinates
  ""
  ([stoin]
     (let [newsto (fs/replace-type stoin "-new.sto")
           fna (sto->fna stoin (fs/replace-type stoin ".fna"))

           ids (fs/replace-type stoin "-id.txt")
           idseq (map #(let [[nm s e st] (-> % first entry-parts flatten)]
                         [nm (str "/" s "-" e "/" st) (second %)])
                      (read-seqs stoin :info :both))
           _ (io/with-out-writer ids (doseq [[id] idseq] (println id)))

           misc (str "-seqidlist " ids " -perc_identity 100.0")
           blastout (blast :blastn fna :strand :both :misc misc)

           id-info (filter-base-hits-to-most-likely blastout)

           [good bad] (entry-aln-sqs id-info idseq)

           bad-file (fs/replace-type stoin "-bad.txt")
           bad (keep (fn[[nm-coord aln-sq]]
                       (let [nm (first (str/split #"/" 2 nm-coord))
                             fna (fs/join default-genome-fasta-dir
                                          (str nm ".fna"))
                             degapped (degap-seqs aln-sq)
                             actual (if (not (fs/exists? fna))
                                      "Non current NC_* sequences unavailable"
                                      (-> nm-coord gen-name-seq second))
                             dup? (= (norm-elements degapped) actual)]
                         [nm-coord aln-sq degapped actual dup?]))
                     bad)

           diff-file (fs/replace-type stoin "-diffs.txt")
           diffs (keep (fn[[nm-coord [ncoord ocoord aln-sq]]]
                         (when (not= ncoord ocoord)
                           [nm-coord ocoord]))
                       good)
           sto-n-orig-line (take 3 (io/read-lines stoin))]

       (io/with-out-writer newsto
         (println (first sto-n-orig-line))
         (println (second sto-n-orig-line))
         (println (third sto-n-orig-line) "\n")
         (doseq [[idi [_ _ sq]] good]
           (cl-format true "~A~40T~A~%" idi sq))
         (doseq [gcline (filter #(.startsWith % "#")
                                (first (sto-GC-and-seq-lines stoin)))]
           (let [[gc kind v] (str/split #" +" 3 gcline)]
             (cl-format true "~A~40T~A~%" (str gc " " kind) v)))
         (println "//"))
       (fs/rename stoin (fs/replace-type stoin "-old.sto"))
       (fs/rename newsto stoin)

       (when (seq diffs)
         (io/with-out-writer diff-file
           (doseq [[nm-coord ocoord] diffs]
             (let [nm (first (str/split #"/" nm-coord))]
               (println (str nm ocoord) " --> " nm-coord)))))

       (when (seq bad)
         (io/with-out-writer bad-file
           (doseq [[nm-coord gapped degapped actual dup?] bad]
             (println nm-coord (if dup? "is DUP!" "") "\n"
                      gapped "\n" degapped "\n"
                      actual "\n"))))

       (when (fs/exists? ids) (io/delete-file ids))
       (when (fs/exists? blastout) (io/delete-file blastout))
       (when (fs/exists? fna) (io/delete-file fna))
       [stoin diff-file bad-file]))

  ([sto1 sto2 & stos]
     (map correct-sto-coordinates (->> stos (cons sto2) (cons sto1)))))


(defn correct-dir-stos
  ""
  [stodir & {:keys [dirdir] :or {dirdir false}}]
  (if (not dirdir)
    (fs/dodir stodir #(fs/directory-files % ".sto") correct-sto-coordinates)
    (fs/dodir stodir #(fs/directory-files % "") correct-dir-stos)))


;;; CD Hit Redundancy removal.
;;;
(defn get-cd-hit-infile [x]
  (if (coll? x)
    (nms-sqs->fasta-file
     (map (fn[[nm rng sq]] [(str nm ":" rng) sq]) x)
     (fs/tempfile "cdhit-" ".fna"))
    (fs/fullpath x)))

(defn cd-hit-est [input outfile & {c :c n :n :or {c 0.90 n 8}}]
  (let [infile (get-cd-hit-infile input)
        outfile (fs/fullpath outfile)
        cdhit-path (get-tool-path :cdhit)
        cdhitcmd (str cdhit-path "cd-hit-est")
        cmdargs ["-i" infile "-o" outfile
                 "-c" (str c) "-n" (str n)]]
    (assert-tools-exist [cdhitcmd])
    (catch-all (runx cdhitcmd cmdargs))
    (when (fs/exists? (str outfile ".clstr"))
      (fs/rm (str outfile ".clstr")))
    (when (fs/exists? (str outfile ".bak.clstr"))
      (fs/rm (str outfile ".bak.clstr")))
    (when (coll? input) (fs/rm infile))
    outfile))




;;; ----------------------------------------------------------------------
;;;
;;; OK, here we have an effort to map EMBL seq/ids to refseq ids.  In
;;; particular, we turn embl tagged stos into nc tagged stos.  So, for
;;; example, RFAM stos can become "refseqed".


(defn get-embl-blast-candidates
  "Takes an ncbi blast csv output file result from blasting an embl
   entry input and returns an embl-id to ncbi-ids map.  The csv has
   fields:

     qseqid qstart qend evalue sseqid start send

   The return map has entries: [embl-id [vec-of-best-ncbi-entries]].
  "
  [blast-output-file]
  (reduce
   (fn[m x]
     (let [embl-entry (first x)
           nc-nm (->> (nth x 4) (str/split #"\|") last
                      (re-find #"^[A-Za-z_0-9]+"))
           ev (Double. (nth x 3))]
       (if (and (= (second x) "1")
                ;; somewhat arbitrary, but real examples should
                ;; have very low evalues
                (< ev 1e-05))
         (let [s (Integer. (nth x 5))
               e (Integer. (nth x 6))
               [s e sd] (if (> s e) [e s -1] [s e 1])
               entry (str nc-nm "/" s "-" e "/" sd)]
           (assoc m embl-entry (conj (get m embl-entry []) entry)))
         m)))
   {} (butlast (csv/parse-csv (slurp blast-output-file)))))


(defn embl-sto->nc-sto
  "Convert EMBL entry names to NCBI entry names while maintaining
   alignment information.  STOIN is a sto file with entries named with
   EMBL accesions.  Result is an output sto file with the alignments
   of STOIN with NCBI NC accensions replacing the EMBL names.  Assumes
   a nucleotide alphabet (no AA support yet).  Uses blast with 100%
   identity match over full database and selects the hits with exact
   size and content match.  If multiple candidates, selects randomly.
  "
  [stoin & {:keys [blastout]}]
  (let [fna (fs/replace-type stoin ".fna")
        ncsto (fs/replace-type stoin "-NC.sto")
        misc (str "-perc_identity 100.0")
        _ (sto->fna stoin fna)
        blastout (or blastout (blast :blastn fna :strand :both :misc misc))
        embl-aln-pairs (->> stoin (#(read-seqs % :info :both)) (into {}))
        nc-aln-pairs (reduce
                      (fn[V [id v]]
                        (conj V [(rand-nth v) (embl-aln-pairs id)]))
                      [] (get-embl-blast-candidates blastout))
        [sqs-gcs hdgf-lines] (sto-GC-and-seq-lines stoin)
        gc-lines (filter #(.startsWith % "#") sqs-gcs)]
    (io/with-out-writer ncsto
      (println (first hdgf-lines))
      (println (second hdgf-lines) "\n")
      (doseq [l (drop 2 hdgf-lines)] (println l))
      (println "\n")
      (doseq [[nc sq] nc-aln-pairs] (cl-format true "~A~40T~A~%" nc sq))
      (doseq [gcline gc-lines]
        (let [[gc kind v] (str/split #"\s+" 3 gcline)]
          (cl-format true "~A~40T~A~%" (str gc " " kind) v)))
      (println "//"))
    ncsto))


(defn embl-to-nc
  "Convert EMBL entry names to NCBI entry names while maintaining
   sequence information.  FIN is either a sto or fna file with entries
   named with EMBL accesions.  Result is a corresponding output sto or
   fna file with the same sequence (and alignment if sto) with NCBI NC
   accensions replacing the EMBL names.  Assumes a nucleotide
   alphabet (no AA support yet).  Uses blast with 100% identity match
   over full database and selects the hits with exact size and content
   match.  If multiple candidates, selects randomly.
  "
  [fin]
  {:pre [(in (fs/ftype fin) ["sto" "fna"])]}
  (if (= (fs/ftype fin) "sto")
    (embl-sto->nc-sto fin)
    (let [fna (fs/fullpath fin)
          ncfna (fs/replace-type fna "-NC.fna")
          misc (str "-perc_identity 100.0")
          blastout (blast :blastn fna :strand :both :misc misc)
          embl-seq-pairs (->> fna (#(read-seqs % :info :both)) (into {}))
          nc-seq-pairs (reduce
                        (fn[V [id v]]
                          (conj V [(rand-nth v) (embl-seq-pairs id)]))
                        [] (get-embl-blast-candidates blastout))]
      (io/with-out-writer ncfna
        (doseq [[nc sq] nc-seq-pairs] (cl-format true ">~A~%~A~%" nc sq)))
      ncfna)))




;;; ----------------------------------------------------------------------
;;;

(defn entry-file-intersect
  "Intersect the entry sets in the files f1 and f2 or f1, f2 and
   remaining fs.  Full is true or false.  If true, entries must fully
   match (name, start, end and strand) to be included in intersection.
   If false, only the names are used to match.  Returns the resulting
   set of entries sans sequences.  See entry-file-intersect-with-seqs
   to include sequences as well.  File types may be fasta (fna,
   hitfna), sto, aln, and gaisr csv.
  "
  ([full f1 f2]
     (let [f1-ncs (read-seqs f1 :info :names)
           f1-ncs (if full
                    (set f1-ncs)
                    (->> f1-ncs
                         (map #(first (entry-parts %)))
                         set))
           f2-ncs (read-seqs f2 :info :names)
           f2-ncs (if full
                    (set f2-ncs)
                    (->> f2-ncs
                         (map #(first (entry-parts %)))
                         set))]
       (set/intersection f1-ncs f2-ncs)))
  ([full f1 f2 & fs]
     (reduce (fn[S f]
               (let [f-ncs (read-seqs f :info :names)
                     f-ncs (if full
                             (set f-ncs)
                             (->> f-ncs
                                  (map #(first (entry-parts %)))
                                  set))]
                 (set/intersection S f-ncs)))
             #{} (conj fs f1 f2))))


(defn entry-file-intersect-with-seqs
  "Like entry-file-intersect, but returns corresponding sqs for the
   result entry set as well.  Uses the full entry (full is true for
   entry-file-intersect).  Intersect the entries in files f1 f2, and
   if non empty those in fs as well.  Take the result set and create
   the set of pairs [entry sq] for each entry in the set.  Returns the
   resulting sequence of pairs.
  "
  [f1 f2 & fs]
  (let [efi-set (apply entry-file-intersect true f1 f2 fs)
        sqf (->> [f1 f2]
                 (filter #(not (in (fs/ftype %) ["ent" "csv"])))
                 first)
        id-sq-map (when sqf (->> sqf (#(read-seqs % :info :both)) (into {})))]
    (reduce (fn[V entry]
              (conj V [entry (if id-sq-map
                               (id-sq-map entry)
                               (gen-name-seq entry))]))
            [] efi-set)))




;;; ----------------------------------------------------------------------
;;; CMFinder tools and operations

(defn process-blast-for-cmfinder
  "A simple reducer for blastn output where the blast run specified
   output format with the options:

   -outfmt '10 qseqid qstart qend evalue sseqid start send'"
  [acc line]
  (let [bs (str/split #"," line)
        qs (take 4 bs)
        ss (drop 4 bs)]
    (conj
     acc
     {:qs (conj (rest qs) (str/replace-re #"^gi\|\d+\|" "" (first qs)))
      :ss (conj (rest ss) (str/replace-re #"^gi\|\d+\|" "" (first ss)))})))




(defn parse-blast
  "Parse a (NCBI) blastn output for use by cmfinder tools.  The blast run
   specified output format with the options:

     -outfmt '10 qseqid qstart qend evalue sseqid start send <& others>'

   IN is the blast output filespec (path or file obj) and OUT is a filespec
   (path or file obj) of where to place the parsed output.  This file
   is then intended to be used as input to the CANDS program of
   cmfinder that generates per motif-loop data files for input to
   cmfinder.
  "
  [in out]
  (do-text-to-text
   [in out]
   (let [bs (str/split #"," $line)
         qseq (take 4 bs)
         sseq (take 3 (drop 4 bs))
         [qid qs qe e] qseq
         [sid ss se] sseq]
     (when (and (not (= qid sid)) (< (Float. e) 0.5))
       (printf "%S :\t %S - %S\t %S\n"
               (str/replace-re #"^gi\|\d+\|" "" qid) qs qe e)
       (printf "%S :\t %S - %S\n\n"
               (str/replace-re #"^gi\|\d+\|" "" sid) ss se))))
  out)




(defnk candf [in out
             :max-cand-seq 40
             :min-len-cand 30
             :max-len-cand 100
             :min-stem-loops 1
             :max-stem-loops 1]

  {:pre [(<= min-stem-loops max-stem-loops)
         (<= max-len-cand 100)
         (<= 1 min-len-cand)
         (<= min-len-cand max-len-cand)]}

  (let [seqin        (fs/fullpath in)
        stem-loops   (if (= max-stem-loops min-stem-loops)
                       (str min-stem-loops)
                       (str min-stem-loops "." max-stem-loops))
        seqout       (or out (str seqin ".h" stem-loops ".cand"))

        cmf-path (get-tool-path :cmfinder)
        candf (str cmf-path "candf")

        args (map str ["-c" max-cand-seq "-o" seqout
                       "-M" max-len-cand "-m" min-len-cand
                       "-s" min-stem-loops "-S" max-stem-loops
                       seqin])]
    (assert-tools-exist [candf])
    (runx candf args)))




(defnk cands [in candf-output
             :max-out-motifs 3
             :expected-motif-freq 0.80
             :blast-seq-match nil]
  (let [seqin (fs/fullpath in)
        cmf-path (get-tool-path :cmfinder)
        cands (str cmf-path "cands")
        args (keep #(when % (str %))
                   (flatten ["-n" max-out-motifs "-f" expected-motif-freq
                             (when blast-seq-match ["-m" blast-seq-match])
                             seqin candf-output]))]
    (assert-tools-exist [cands])
    (runx cands args)))




(defn canda [seqin cands-file canda-ofile]
  (let [cmf-path (get-tool-path :cmfinder)
        canda (str cmf-path "canda")]
    (assert-tools-exist [canda])
    (runx canda
          (fs/fullpath cands-file)
          (fs/fullpath seqin)
          (fs/fullpath canda-ofile))))




(defn cmfinder
  [seqin canda-file motif-out cm-out
   & {:keys [initial-cm candf-output stdout]}]

  (let [cmf-stdout   (and stdout (fs/fullpath stdout))
        initial-cm   (and initial-cm (fs/fullpath initial-cm))
        candf-output (and candf-output (fs/fullpath candf-output))

        cmf-path (get-tool-path :cmfinder)
        cmfinder (str cmf-path "cmfinder")

        args (keep #(when % %)
                   (flatten ["-o" (fs/fullpath motif-out)
                             "-a" (fs/fullpath canda-file)
                             (when initial-cm ["-i" initial-cm])
                             (when candf-output ["-c" candf-output])
                             (fs/fullpath seqin)
                             (fs/fullpath cm-out)
                             (when cmf-stdout [:> cmf-stdout])]))]
    (assert-tools-exist [cmfinder])
    (runx cmfinder args)))




(defn cmfinder*
  "Perform a full run of cmfinder sub pipeline.  Basically:

   (-> in-fna (candf options)
       (#(if blastpgm) (blast options))
       (cands options)
       (map #(-> % (canda in-fna) (cmfinder in-fna))))

   Where options are the given set of optional keys and values (with
   noted defaults).
  "
  [in & {:keys [blastpgm blastdb initial-cm
                strand word-size max-cand-seq
                min-len-cand max-len-cand
                max-out-motifs
                expected-motif-freq
                min-stem-loops max-stem-loops
                del-files]
         :or {strand :plus
              word-size 8
              max-cand-seq 40
              min-len-cand 30
              max-len-cand 100
              max-out-motifs 3
              expected-motif-freq 0.80
              min-stem-loops 1
              max-stem-loops 1
              del-files true}}]

  (let [seqin        (fs/fullpath in)
        stem-loops   (if (= max-stem-loops min-stem-loops)
                        (str min-stem-loops)
                        (str min-stem-loops "." max-stem-loops))
        seqout-candf (str seqin ".h" stem-loops ".cand")
        seq-match    (str seqin ".match")]

    ;; CANDF --> produce stem loop candidates (all in seqout-candf)
    (candf seqin seqout-candf
           :max-cand-seq max-cand-seq
           :min-len-cand min-len-cand :max-len-cand max-len-cand
           :min-stem-loops min-stem-loops :max-stem-loops max-stem-loops)

    ;; When blasting, produce blast match file SEQ-MATCH, for input to cands
    (when blastpgm
      (blast  blastpgm seqin (str seqin ".blast")
              :blastdb blastdb
              :strand strand :word-size word-size
              :fn parse-blast :fnargs seq-match))

    ;; CANDS --> produce MAX-OUT-MOTIF candidate files for CANDA/CMFINDER
    (cands seqin seqout-candf
           :max-out-motifs max-out-motifs
           :expected-motif-freq expected-motif-freq
           :blast-seq-match (when blastpgm seq-match))

    ;; Now, for each motif candidate MC, canda MC ... | cmfinder ...
    ;; which yields an alignment and motif (hits) set of output.
    (loop [i (int 0)
           motif-stos []]
      (if (= i max-out-motifs)
        motif-stos
        (let [dot-i (str "." (+ 1 i))
              suffix (str stem-loops dot-i)
              cands-file (str seqout-candf dot-i)]
          (when (not (fs/empty? cands-file))
            (let [canda-ofile (str seqin ".align-sto.h" suffix)
                  motif-ofile (str seqin ".motif-sto.h" suffix)
                  cm-file     (str seqin ".cm.h"    suffix)
                  cmf-stdout  (str seqin ".h" stem-loops ".out" dot-i)]
              (canda seqin cands-file canda-ofile)
              (cmfinder seqin canda-ofile
                        motif-ofile cm-file
                        :initial-cm initial-cm
                        :stdout cmf-stdout)
              (when del-files
                (fs/rm cands-file)
                (fs/rm canda-ofile)
                (fs/rm cm-file)
                (fs/rm cmf-stdout))
              (recur (unchecked-inc i)
                     (conj motif-stos motif-ofile)))))))
    ))




(defn cmalign
  "Run an automated parallelized cmalign on originating covariant
   model CM (presumably obtained by some prior cmbuild/calibrate) and
   fasta file SEQFILE.  Place new output sto in OUTFILE.  Return
   outfile.  OPTS is any set of k/v pairs (as a vector) that are legal
   for cmalign.

   This is mostly (only?) used when building the next version of an
   input sto that has gone through a complete cmbuild through cmsearch
   through FFP and/or Scribl process yielding a new set of candidate
   targets that have been placed in seqfile.
  "
  [cm seqfile outfile
   & {opts :opts par :par :or {opts ["-q" "-1"] par 3}}]
  (let [infernal-path (get-tool-path :infernal)
        cmaligncmd (str infernal-path "cmalign")
        cmfile (fs/fullpath cm)
        seqfile (fs/fullpath seqfile)
        outfile (fs/fullpath outfile)
        mpirun "mpirun"
        cmdargs (conj (vec (concat ["-np" (str par)
                                    cmaligncmd "--mpi" "-o" outfile]
                                   opts))
                      cmfile seqfile)]
    (if (fs/empty? seqfile)
      (raise :type :empty-seq-file :file seqfile)
      (do (assert-tools-exist [cmaligncmd])
          (runx mpirun cmdargs)))
    outfile))
