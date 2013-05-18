;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                             S E Q - U T I L S                            ;;
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

(ns edu.bc.bio.seq-utils

  "Intended to include/contain various ancillary operations on
   sequence data. Refactoring will sooner than later make this
   obsolete.  Note it now contains what once was in seq-utils2 and
   that should be factored as well."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        [clojure.pprint
         :only [cl-format]]
        ))


;;; ----------------------------------------------------------------------
;;;
;;; IUPAC (International Union of Pure and Applied Chemistry) codes
;;; for nucleotides and various "groups" of nucleotides.
;;;

(def +IUPAC+
     ^{:doc "These are the codes for \"bases\" used in alignments"}
     {\A "Adenine"
      \C "Cytosine"
      \G "Guanine"
      \T "Thymine"
      \U "Uracil"
      \R "AG"
      \Y "CUT"
      \S "GC"
      \W "AUT"
      \K "GUT"
      \M "AC"
      \B "CGUT"
      \D "AGUT"
      \H "ACUT"
      \V "ACG"
      \N "any"
      \. "gap"
      \- "gap"})

(def +NONSTD-RNA+
     ^{:doc "These are the codes for all non standard \"bases\" used
             in alignments"}

     {\T "Thymine"
      \R "AG"       ; purines
      \Y "CUT"      ; pyrimidines
      \S "GC"
      \W "AUT"
      \K "GUT"
      \M "AC"
      \B "CGUT"
      \D "AGUT"
      \H "ACUT"
      \V "ACG"
      \N "any"
      \. "gap"
      \- "gap"})

(def +RY-XLATE+
     ^{:doc "Map for translating DNA and RNA seqs to purine and
             pyrimidine seqs.  Assumes normed-elements (for which
             see)"}
     {\A \R, \G \R, \C \Y, \U \Y, \T \Y})


;;; ----------------------------------------------------------------------
;;;
;;; Various sequence (element) transformers and translators.  Mostly
;;; for nucleotide sequences, but could (should?) be generalized a bit
;;; for AA seqs as well.
;;;


(defn norm-elements
  "\"Normalize\" elements in sequences by ensuring each character is
   mapped to its uppercase variant."
  [seqs]
  (if (string? seqs)
    (str/upper-case seqs)
    (map str/upper-case seqs)))


(defn gap-percent
  "Return the percentage of sq (typically an aligned variant of a
   genome sequence) that is comprised of gap characters.  For the
   single parameter case, gaps are taken as \\. and \\-.  For the
   gap-chars case, gap characters are the elements (characters) in
   gap-chars (a seqable collection)."
   ([sq]
      (let [[_ ps] (freqs-probs 1 sq)
            gp (+ (get ps \- 0) (get ps \. 0))]
        gp))
   ([sq gap-chars]
      (let [[_ ps] (freqs-probs 1 sq)
            gp (sum (filter #(get ps % 0) gap-chars))]
        gp)))

(defn filter-pgap
  "Filter the sequence set seqs, by returning only those that have
   less than a pgap percentage of gap characters. Gap chars are either
   the defaults \\. \\- or those in gap-chars (a seqable collection)."
  ([seqs pgap]
     (filter #(< (gap-percent %) pgap) seqs))
  ([seqs pgap gap-chars]
     (filter #(< (gap-percent % gap-chars) pgap) seqs)))

(defn degap-seqs
  "Remove gap characters from a sequence or set of sequences.  These
   would be from an alignment set.  It is not clear how / where useful
   this is/would be as it destroys the alignment spacing.  Other than
   a single sequence situation, use degap-tuples instead!!

   Gap chars are either the defaults \\. \\- or those in gap-chars (a
   seqable collection).
  "
  ([seqs]
     (if (coll? seqs)
       (map #(str/replace-re #"[-.]+" "" %) seqs)
       (str/replace-re #"[-.]+" "" seqs)))
  ([seqs gap-chars]
     (if (coll? seqs)
       (map (fn[sq] (apply str (filter (fn[c] #(not (in c gap-chars))) sq)))
            seqs)
       (apply str (filter (fn[c] #(not (in c gap-chars))) seqs)))))

(defn gaps?
  "Return whether K, a char, string, or coll, contains a \"gap
   character\".  For the single parameter case gap chars are taken as,
   a \\. or \\-.  For the gap-chars case gap characters are the
   elements (characters) in gap-chars (a seqable collection).
  "
  ([k]
     (or
      (and (char? k) (in k [\. \-]))
      (and (string? k) (re-find #"(\.|-)" k))
      (and (coll? k) (or (in \- k) (in \. k)))))
  ([k gap-chars]
     (or
      (and (char? k) (in k gap-chars))
      (and (string? k) (some #(in % gap-chars) k))
      (and (coll? k) (some #(in % gap-chars) k)))))


(defn gap-count
  "Return the number of \"gap characters\" in sq.  For the single
   parameter case gap chars are taken as, a \\. or \\-.  For the
   gap-chars case gap characters are the elements (characters) in
   gap-chars (a seqable collection).
  "
  ([sq]
     (gap-count sq [\. \-]))
  ([sq gap-chars]
     (reduce (fn[cnt ch] (+ cnt (if (gaps? ch gap-chars) 1 0))) 0 sq)))




(defn reverse-compliment
  "Reverse compliment of DNA/RNA sequence SQ. Return the sequence that
   would pair with (reverse SQ).

   Ex: (reverse-compliment \"AAGGAAUUCC\") => \"GGAAUUCCUU\"
  "
  [sq]
  (let [ucp (first (str/codepoints "U"))
        m (if (neg? (.indexOf sq ucp))
            {\A \T, \T \A, \G \C, \C \G}
            {\A \U, \U \A, \G \C, \C \G})
        ;; Add in the 'odd balls'
        m (merge m {\R \Y, \Y \R,
                    \K \M, \M \K,
                    \B \V, \V \B,
                    \D \H, \H \D,
                    \S \S, \W \W, \N \N})]
    (str/join "" (reverse (map m sq)))))


(defn seqXlate
  "Translate sqs according to the translation map xmap.  SQS is either
   a single sequence (string) or a collection of such.  It cannot be a
   file.  Each sequence has its elements transformed according xmap
   and each is returned as a string.
  "
  [sqs & {:keys [xmap] :or {xmap +RY-XLATE+}}]
  (let [xlate (fn [sq]
                (apply str (map #(let [x (xmap %)] (if x x %))
                                (degap-seqs sq))))]
    (if (string? sqs)
      (xlate sqs)
      (map xlate sqs))))




;;; ----------------------------------------------------------------------
;;;
;;; Sequence shuffling with nucleotide/aa frequency distributions
;;; preserved.  Currently only provides support for dinucleotide
;;; distributions.  Which is basically bogus, but the most typical use
;;; case??


(defn nt-cntn
  "Frequencies of N-(nucleotides/amino-acids) of sequence SEQUENCE,
   i.e., for N=2 and bases freq of dinucleotides.  Effectively a
   rename of utils:freqn.
   "
  [n sequence]
  (freqn n sequence))


(def test-nts
     (nt-cntn 2 (str/join "" (repeat 1 "acagtcaacctggagcctggt"))))

(def nts-map
     (reduce #(assoc %1 (subs (first %2) 0 1)
                     (conj (get %1 (subs (first %2) 0 1) []) %2))
             {} test-nts))

(def legal-ends
     (keep #(when (= (subs (first %1) 1 2) "t") %1) test-nts))


(def limit-counter (atom 50000))
(def print-limit-info (atom nil))

(defn check-limit [start len ss newseq nts-map]
  (swap! limit-counter dec)
  (when (= limit-counter 0)
    (when print-limit-info
      (prn start len ss newseq nts-map))
    (raise :type :explode)))


(def shuffles nil)

(defn dint-shuffle [start len nts-map legal-ends ss newseq & base-starters]
  (cond
   (and (= len 0) (some #(= (last newseq) (first %1)) legal-ends))
   (str ss (subs (last newseq) 1 2))

   (= len 0) :fail

   (nil? (nts-map start)) :fail

   :otherwise
   (loop [starters (or (first base-starters) (nts-map start))
          [pick cnt] (first starters)]
     (let [newstarters (drop 1 starters)
           mod-edges (keep  #(if (not= pick (first %1)) %1
                               (let [k (first %1) cnt (dec (second %1))]
                                 (when (> cnt 0) [k cnt])))
                            (nts-map start))
           new-map (if (not-empty mod-edges)
                     (assoc nts-map start
                            (concat (drop 1 mod-edges)
                                    (take 1 mod-edges)))
                     (dissoc nts-map start))
           result (dint-shuffle (subs pick 1 2) (dec len) new-map
                                legal-ends (str ss start)
                                (conj newseq pick))]
       (if (and result (not= result :fail)
                (not (@shuffles result)))
         result
         (do
           (check-limit start len ss newseq nts-map)
           (if (not-empty newstarters)
             (recur newstarters
                    (first newstarters))
             :fail)))))))


(defn dint-seq-shuffle [limit s]
  (let [len (count s)
        start (subs s 0 1)
        end (subs s (dec len) len)
        di-nts (nt-cntn 2 s)
        nts-map (reduce #(assoc %1 (subs (first %2) 0 1)
                                (conj (get %1 (subs (first %2) 0 1) []) %2))
                        {} di-nts)
        legal-ends (keep #(when (= (subs (first %1) 1 2) end) %1) di-nts)]
    (binding [shuffles (atom {})
              limit-counter limit-counter]
      (handler-case
        (loop [_ limit
               base-starters (nts-map start)]
          (let [s (dint-shuffle
                   start (dec len) nts-map legal-ends "" [] base-starters)]
            (when (not= :fail s)
              (swap! shuffles
                     assoc s (inc (get shuffles s 0))))
            (when (> _ 1)
              (recur (dec _)
                     (concat (drop 1 base-starters)
                             (take 1 base-starters))))))
        ([:type :explode] e
           (println :limit-reached)))
        @shuffles)))




;;; Simple Monte Carlo subseq frequency approximation within total seqs
;;;

(defn alphabet [type]
  {:pre [(in type [:dna :rna :amino])]}
  (case type
        :dna ["A" "T" "G" "C"]
        :rna ["A" "U" "G" "C"]
        :amino ["A" "R" "N" "D" "C" "Q" "E" "G" "H" "I"
                "L" "K" "M" "F" "P" "S" "T" "W" "Y" "V"]))


(defn gen-random-bioseq  [type size & seed]
  {:pre [(> size 0)]}
  (let [alphabet (alphabet type)
        seed (first seed)
        seed (if seed (vec (str/upper-case seed)) [])
        size (- size (count seed))]
    (loop [s seed
           cnt 0]
      (if (< cnt size)
        (recur (conj s (rand-nth alphabet))
               (inc cnt))
        [(str/join "" s) s]))))


(defn monte-carlo-subseq-cnts-body [type sample-size seqsize sseq len seed]
  (loop [ss sample-size
         cnts {}
         sseqcnts {}]
    (if (= ss 0)
      [cnts sseqcnts]
      (let [s (first (gen-random-bioseq type seqsize (when seed sseq)))
            scnts (nt-cntn len s)]
        (recur (dec ss)
               (reduce (fn[m [k v]]
                         (assoc m k (+ (get m k 0) v)))
                       cnts scnts)
               (assoc sseqcnts (scnts sseq) (inc (get scnts sseq 0))))))))


(defnk monte-carlo-subseq-cnts
  [:type :dna :size 1024 :sseq "" :sample-size 10000 :seed nil]
  {:pre [(> size 0)
         (let [abet (alphabet type)]
           (every? #(in (str %1) abet) (str/upper-case sseq)))]}
  (let [len (count sseq)
        sseq (str/upper-case sseq)]
    (monte-carlo-subseq-cnts-body type sample-size size sseq len seed)))


(defn reduce-seqcnt-results [cntmaps]
  (loop [cntmaps cntmaps
         total-cnts {}
         sseqcnts {}]
    (if (empty? cntmaps)
      [total-cnts sseqcnts]
      (let [[tcnts scnts] (first cntmaps)]
        (recur (drop 1 cntmaps)
               (reduce (fn[m [k v]]
                         (assoc m k (+ (get m k 0) v)))
                       total-cnts tcnts)
               (reduce (fn[m [k v]]
                         (assoc m k (+ (get m k 0) v)))
                       sseqcnts scnts))))))

(defnk pmonte-carlo-subseq-cnts
  [:type :dna :size 1024 :sseq "" :sample-size 10000 :seed nil :par 10]
  {:pre [(> size 0)
         (and (integer? par) (<= 1 par 12))
         (let [abet (alphabet type)]
           (every? #(in (str %1) abet) (str/upper-case sseq)))]}
  (let [len (count sseq)
        sseq (str/upper-case sseq)
        psz (floor (/ sample-size par))
        [q r] (div sample-size psz)
        ssizes (repeat q psz)
        ssizes (conj (drop 1 ssizes) (+ (first ssizes) r))]
    (reduce-seqcnt-results
     (pmap #(monte-carlo-subseq-cnts-body type %1 size sseq len seed)
           ssizes))))


(defnk monte-carlo-subseq-freq
  [:type :dna :size 1024 :sseq "" :sample-size 10000 :seed nil :par 10]
  (let [f (if (and par (integer? par) (> par 1))
            pmonte-carlo-subseq-cnts
            monte-carlo-subseq-cnts)
        [cnts sseqcnts] (f :type type :size size :sseq sseq
                           :sample-size sample-size :seed seed :par par)
        sseq-total (cnts (str/upper-case sseq))
        total (reduce (fn[cnt [k v]] (+ cnt v)) 0 cnts)]
    [sseq-total total (double (/ sseq-total total)) sseqcnts]))



