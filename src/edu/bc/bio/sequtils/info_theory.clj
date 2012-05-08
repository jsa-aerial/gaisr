;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;               S E Q U T I L S . I N F O - T H E O R Y                    ;;
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

(ns edu.bc.bio.sequtils.info-theory

  "Various information theory computations, calculations and results
   as applied to bio sequences and alignments thereof.  Includes
   entropy, joint entropy, conditional entropy, mutual information,
   conditional mi, et. al."

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
        edu.bc.bio.sequtils.files
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        ))




(defn norm-elements
  "\"Normalize\" elements in sequences by ensuring each character is
   mapped to its uppercase variant."
  [seqs]
  (map str/upper-case seqs))


(defn gap-percent
  "Return the percentage of seq (typically an aligned variant of a
   genome sequence) that is comprised of gap characters.  See gaps?."
   [seq]
  (let [[_ ps] (freqs-probs 1 seq)
        gp (+ (get ps \- 0) (get ps \. 0))]
    gp))

(defn filter-pgap
  "Filter the sequence set seqs, by returning only those that have less
   than a pgap percentage of gap characters."
  [seqs pgap]
  (filter #(< (gap-percent %) pgap) seqs))

(defn gaps?
  "Return whether K, a char, string, or coll, contains a \"gap
   character\", i.e., a . or -.

   Note: bug: may need a more complete set of what these are or
   perhaps pass the set in with a default of [. -].
  "
  [k]
  (or
   (and (char? k) (in k [\. \-]))
   (and (string? k) (re-find #"(\.|-)" k))
   (and (coll? k) (or (in \- k) (in \. k)))))

(defn degap-freqs
  "Take a frequency map, and remove all elements whose keys contain
   gap characters ('-', '.')
  "
  [freq-map]
  (reduce (fn[m [k v]]
            (if (gaps? k) m (assoc m k v)))
          {} freq-map))

(defn degap-seqs
  "Remove gap characters from a sequence or set of sequences.  These
   would be from an alignment set.  It is not clear how / where useful
   this is/would be as it destroys the alignment spacing
  "
  [seqs]
  (if (coll? seqs)
    (map #(str/replace-re #"[-.]+" "" %) seqs)
    (str/replace-re #"[-.]+" "" seqs)))




(defn seqs-shannon-entropy
  "Returns the Shannon entropy of the set of sequences in SEQSET, a
   collection of sequences or a string denoting a legal format
   sequence file (see read-seqs).  Returns a pair [ses total-ses],
   where

   ses is (map shannon-entropy seqset) and total-ses is the total over
   all of seqset.
  "
  [seqset]
  {:pre [(gaisr-seq-set? seqset)]}
  (let [ses (if (coll? seqset)
              (map shannon-entropy seqset)
              (map shannon-entropy (read-seqs seqset)))]
    [ses (/ (sum ses) (count ses))]))


(defn seqs-freqs-probs
  "Return sequence frequencies and probabilities over the set of
   sequences in SEQSET, a collection of sequences or a string denoting
   a legal format sequence file (see read-seqs).  FSPS-FN is a
   function taking a combination count or window width of N and a
   sequence collection (here, seqset) and optional par parameter.
   Applies fsps-fn to n and seqset.

   For large (count seqset) with expensive fss-fn, use par to
   parallelize computation over par chunks.

   Return [ccfsps allfs allps tcount], where

   ccfs&ps is a seq of triples [fs ps cnt], for each C in colls,
   allfs is the map of freqs over all Cs,
   allps is the map of probs over all Cs and
   tcount is the total items over coll.

   NOTE: item keys are \"stringified\", i.e., if k is a key from a map
   produced by fsps-fn, ensures that all returned maps use (apply str
   k) for all keys.
  "
  [n seqset &
   {fsps-fn :fsps-fn nogaps :nogaps par :par
    :or {fsps-fn cc-freqs-probs nogaps true par 1}}]
  {:pre [(gaisr-seq-set? seqset)]}
  (let [seqset (if (coll? seqset) seqset (read-seqs seqset))
        rfn #(reduce (fn[m [k v]] (assoc m (apply str (ensure-vec k)) v)) {} %)
        str-em (fn[[fs ps cnt]] [(rfn fs) (rfn ps) cnt])

        [ccfsps tfs tps tcnt] (if (= par 1) ; Accomodate fsps-fns w/o par arg!
                                (fsps-fn n seqset)
                                (fsps-fn n seqset :par par))

        ccfsps (pxmap str-em par ccfsps)]
    [ccfsps (rfn tfs) (rfn tps) tcnt]))


(defn bg-freqs-probs "NOT YET IMPLEMENTED" [] :NYI)




(defn aln-freqs-probs
  "Take the sequences in SEQSET, a collection of sequences or a string
   denoting a legal format sequence file (see read-seqs), treat as a
   matrix M encoding an alignment (and so all seqs in set must have
   equal length).  If NORM, normalize all base characters to uppercase
   in M.  PGAP is the percent of gaps cutoff for columns.  Filter M,
   by removing all columns Ci, (> (gap-percent Ci) pgap) to get M'.

   Uses FSPS-FN, a function taking a combination count or window width
   of N and a sequence collection, to compute the frequencies and
   probabilities of columns Ci in M' by means of seqs-freqs-probs.
   Let Ci-fs-ps be the results for Ci.  Typically such a result would
   be a triple [fs ps cnt], where fs and ps are maps of frequencies
   and corresponding probabilities keyed by the n-tuple of bases (and
   possibly gaps) underlying fsps-fn (see, for example,
   cc-freqs-probs).

   If NOGAPS is true, remove all items with gaps from maps and
   recompute new probabilities for resulting reduced frequencies sets.

   For large (count seqset) with expensive fss-fn, use par to
   parallelize computation over par chunks.

   Returns [ccfs&ps allfs allps tcount], where

   ccfs&ps is the seq Ci-fs-ps calculated from M'
   allfs is the map of freqs obtained by reducing over all Ci-fs maps
   allps is the map of probs obtained by reducing over all Ci-ps maps
   tcount is the total item (obtained n-tuples) count.
  "
  [n seqset & {fsps-fn :fsps-fn par :par nogaps :nogaps norm :norm pgap :pgap
               :or {fsps-fn cc-combins-freqs-probs par 1
                    nogaps true norm true pgap 0.70}}]

  {:pre [(gaisr-seq-set? seqset)]}

  (let [aln (if (coll? seqset)
              (transpose seqset)
              (read-aln-seqs seqset :cols true))
        alnx (filter-pgap aln pgap)
        alny (if norm (norm-elements alnx) alnx)

        [colsfsps fs ps cnt :as all]
        (seqs-freqs-probs n alny :fsps-fn fsps-fn :par par)

        degap degap-freqs
        reprob (fn[sz fs]
                 (reduce (fn[m [k v]] (assoc m k (double (/ v sz))))
                         {} fs))
        adjust-one (fn[[fs ps cnt]]
                     (let [fs (degap fs)
                           sz (sum fs)
                           ps (reprob sz fs)]
                       [fs ps sz]))]
    (if (not nogaps)
      all
      (let [colsfsps (pxmap adjust-one par colsfsps)
            fs (degap fs)
            sz (sum fs)
            ps (reprob sz fs)]
        [colsfsps fs ps sz]))))


(defn aln-entropy
  "Compute the entropy of each column Ci of an alignment given in SEQSET,
   a gaisr-seq-set.  Entropy is based on the freqs and probs of
   elements of Ci taken n at a time.  ARGS are any keyword arguments
   taken by aln-freqs-probs.

   The manner of taking the elements is determined by the fsps-fn
   argument of aln-freqs-probs.  The default for this,
   cc-combins-freqs-probs is based on the combins function, which
   generates all n-element subsets of Ci.  cc-freqs-probs is based on
   freqn which generates the sliding window of n elements from Ci.

   Returns [cols-entropies total-entropy tcnt], where

   cols-entropies is a seq of entropies for each Ci
   total-entropy is the total entropy over all the columns (from total probs)
   tcnt is the total over all columns of elements counted
  "
  [n seqset & args]
   {:pre [(gaisr-seq-set? seqset)]}
   (let [[col-fs-ps allfs allps tcnt] (apply aln-freqs-probs n seqset args)
         entropy (fn[probs] (- (sum #(* % (log2 %)) probs)))]
     [(map entropy (map #(vals (second %)) col-fs-ps))
      (entropy (vals allps))
      tcnt]))

(defn aln-shannon-entropy
  "Application of aln-entropy with cc-freqs-probs and n=1.  So,
   shannon entropy of each column and totals over all columns.
  "
  [seqset & args]
  (apply aln-entropy 1 seqset :fsps-fn cc-freqs-probs args))

(defn aln-joint-entropy
  "Application of aln-entropy with cc-combins-freqs-probs and n=2.
   So, joint entropy of each column with itself and overall totals.
  "
  [seqset & args]
  (apply aln-entropy 2 seqset :fsps-fn cc-combins-freqs-probs args))




(defn- adjust-seqs-info
  [sqs cols norm pgap]
  (cond
   (and cols norm) (-> sqs transpose (filter-pgap pgap) norm-elements)
   cols (-> sqs transpose (filter-pgap pgap) norm-elements)
   norm (-> sqs (filter-pgap pgap) norm-elements)
   :else (filter-pgap sqs pgap)))

(defn- account-for-symmetry [triples]
  (let [x (coalesce-xy-yx triples
                          (fn [x v] (if (not v) () (cons (second x) v))))]
    (reduce (fn[c [x vs]]
              (reduce (fn[c v]
                        (conj c [x v]))
                      c vs))
            [] x)))

(defn aln-conditional-mutual-information
  ""
  [seqset & {par :par nogaps :nogaps pgap :pgap cols :cols norm :norm
             :or {par 1 nogaps true pgap 0.25 cols true norm true}}]
  {:pre [(gaisr-seq-set? seqset)]}
  (let [aln (if (coll? seqset) seqset (read-aln-seqs seqset))
        aln (adjust-seqs-info aln cols norm pgap)
        aln-map (reduce (fn[m s] (assoc m s true)) {} aln)
        seq-pairs (combins 2 aln)

        triples (pxmap (fn[[x y]]
                         (partition-all
                          2 (interleave
                             (transpose [x y])
                             (map (fn[fs]
                                    (-> fs ((partial freqn 1))
                                        (#(if nogaps (degap-freqs %) %))))
                                  (transpose (keys (dissoc aln-map x y)))))))
                       par
                       seq-pairs)
        triples (if nogaps
                  (pxmap #(filter (fn[[xy _]] (not (gaps? xy))) %) par triples)
                  triples)
        triples (pxmap account-for-symmetry par triples)

        xyz-joint-probs (map #(probs 1 %) triples)
        z-probs (map (fn[tpl] (probs 1 (map second tpl)))
                     triples)
        xz-joint-probs (pxmap (fn[tpl]
                                (probs 1 (map (fn[[xy z]]
                                                [(first xy) z])
                                              tpl)))
                              par
                              triples)
        yz-joint-probs (pxmap (fn[tpl]
                                (probs 1 (map (fn[[xy z]]
                                                [(second xy) z])
                                              tpl)))
                              par
                              triples)
        Ixy|z (pxmap  (fn[xyz-jps z-ps xz-jps yz-jps]
                        (sum (fn[xyz-jp zp xz-jp yz-jp]
                               (* xyz-jp (log2 (/ (* xyz-jp zp)
                                                  (* xz-jp yz-jp)))))
                             :|| xyz-jps z-ps xz-jps yz-jps))
                      par
                      (map vals xyz-joint-probs)
                      (map vals z-probs)
                      (map vals xz-joint-probs)
                      (map vals yz-joint-probs))]
    [Ixy|z xyz-joint-probs xz-joint-probs yz-joint-probs z-probs]))


(defn aln-mutual-information
  ""
  [seqset & {par :par nogaps :nogaps norm :norm pgap :pgap cols :cols
             :or {par 1 nogaps true norm true pgap 0.25 cols false}}]
  {:pre [(gaisr-seq-set? seqset)]}
  (let [probs (fn[fs]
                (let [sz (sum fs)]
                  (reduce (fn[m [k v]] (assoc m k (double (/ v sz))))
                          {} fs)))
        entropy (fn[probs] (- (sum #(* % (log2 %)) (vals probs))))
        aln (if (coll? seqset) seqset (read-aln-seqs seqset))
        aln (adjust-seqs-info aln cols norm pgap)

        seq-entropy-map (into {} (map #(do [% (shannon-entropy
                                               (if nogaps (degap-seqs %) %))])
                                      aln))

        seq-pairs (combins 2 aln)
        pair-freqs (map #(combin-count-reduction % true)
                        (map transpose seq-pairs))
        pair-freqs (if nogaps (map degap-freqs pair-freqs) pair-freqs)

        joint-probs (map probs pair-freqs)
        joint-entropy (map entropy joint-probs)
        mutual-info (map (fn[Hxy [x y]]
                           (let [Hx (seq-entropy-map x)
                                 Hy (seq-entropy-map y)]
                             (+ Hx Hy (- Hxy))))
                         joint-entropy seq-pairs)]
    [mutual-info joint-entropy joint-probs pair-freqs]))


(defn variation-information
  [[seqx seqy] & {nogaps :nogaps :or {nogaps true}}]
  (let [seqs [seqx seqy]
        [[Ixy] [Hxy]] (aln-mutual-information seqs :pgap 1.0 :nogaps nogaps)]
    (cond
     (and (= Hxy 0.0) (not= seqx seqy)) 1
     (= Hxy 0.0) 0
     :else
     (- 1.0 (/ Ixy Hxy)))))
