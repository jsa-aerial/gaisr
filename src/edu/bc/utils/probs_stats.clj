;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                 U T I L S . P R O B S - S T A T S                        ;;
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

(ns edu.bc.utils.probs-stats

  "Various frequency, combinatorial, probability, statistical,
   measure, and metrics for a variety of sequence and string data."

  (:require [clojure.contrib.math :as math]
            [clojure.contrib.combinatorics :as comb]
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        [clojure.contrib.condition
         :only [raise handler-case *condition*
                print-stack-trace stack-trace-info]]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]))


;;; -----------------------------------------------------------------
;;; Various frequency counts, combinations and combinators,
;;; probabilities, similarity coefficients, corresponding difference
;;; fns, ngram operations, edit distances, etc.

(defn letter-pairs [n s]
  (set (map #(apply str %) (partition n 1 s))))

(defn word-letter-pairs [n s]
  (reduce (fn[m s] (set/union m (letter-pairs n s))) #{} (str/split #" " s)))


(defn freqn
  "Frequencies of n-grams in collection COLL treated as a sequence
   Ex: (freqn 2 \"acagtcaacctggagcctggt\")
   =>
   {\"aa\" 1, \"cc\" 2, \"gg\" 2, \"ac\" 2, \"ag\" 2, \"gt\" 2,
   \"tc\" 1, \"ct\" 2, \"tg\" 2, \"ga\" 1, \"gc\" 1, \"ca\" 2}
  "
  [n coll]
  (if (= 1 n)
    (frequencies (seq coll))
    (loop [s (seq coll)
           res (transient {})]
      (let [k (str/join "" (take n s))]
        (if (>= (count s) n)
          (recur (drop 1 s)
                 (assoc! res k (inc (get res k 0))))
          (persistent! res))))))

(defn freqs&probs
  "Generate frequencies of items in coll by means (freq-fn n coll),
   take those and compute the items' corresponding probabilities.
   Return triple: [freqs probs cnt], where freqs and probs are maps
   keyed by item and cnt is the total count of all items.  Lower level
   functional base for more public functions freqs-probs, et.al.
  "
  [n coll freq-fn]
  (let [freqs (freq-fn n coll)
        sz (sum (vals freqs))
        probs (reduce (fn[m [k v]]
                        (assoc m k (float (/ v sz))))
                      {} freqs)]
    [freqs probs sz]))

(defn freqs-probs
  "freqs&probs for n and coll with freq-fn equal freqn"
  [n coll]
  (freqs&probs n coll freqn))

(defn probs
  "Probabilities of items generated from n and coll by freqn"
  [n coll]
  (second (freqs-probs n coll)))

(defn alphabet2
  "Return the alphabet for coll.  If coll is a map, return it's keys.
   If a set, return coll, otherwise return alphabet generated over
   coll by freqn with size n."
  [coll & {n :n :or {n 1}}]
  (cond
   (map? coll) (keys coll)
   (set? coll) coll
   :else
   (keys (freqn n coll))))


(defn cc-freqs
  "Frequencies of items in each collection in the collection of
   collections colls.  Uses freq-fn to generate and count items.  Use
   par to parallelize when (count colls) is large and each (freq-fn
   c), c in colls, is expensive.  Returns a seq of frequency maps, one
   for each collection in colls."
  [n colls freq-fn & {par :par :or {par 1}}]
  (pxmap #(freq-fn n %) par colls))

(defn cc-tfreqs
  "Same as cc-freqs, but also compute total item frequencies reduced
   over all collections in colls.  Returns pair [ccfreqs cctfreqs],
   where ccfreqs is the seq of frequency maps per cc-freqs and
   cctfreqs is the single map for all items.
  "
  [n colls freq-fn & {par :par :or {par 1}}]
  (let [ccfreqs (cc-freqs n colls freq-fn :par par)
        cctfreqs (reduce
                  (fn[m cfreq]
                    (merge-with #(+ %1 %2) m cfreq))
                  {} ccfreqs)]
    [ccfreqs cctfreqs]))

(defn cc-freqs&probs
  "For each C in colls (a collection of collecions), compute the
   frequency and probability map for C and its total item count as
   given by freqs&probs-fn (for example, combins-freqs-probs).  Using
   these maps, compute the total item frequency and probablity maps
   and total count over all collections in colls.  Low level
   functional base for public functions such as cc-freqs-probs,
   et. al.

   For large (count colls) with expensive freqs&probs-fn, use par to
   parallelize computation over par chunks.

   Return [ccfs&ps allfs allps tcount], where

   ccfs&ps is a seq of triples [fs ps cnt], for each C in colls,
   allfs is the map of freqs over all Cs,
   allps is the map of probs over all Cs and
   tcount is the total items over coll.
  "
  [n colls freqs&probs-fn &
  {par :par :or {par 1}}]
  (let [ccfs&ps (pxmap #(freqs&probs-fn n %) par colls)
        allfs (reduce (fn[m xyz] (merge-with + m (first xyz))) {} ccfs&ps)
        sz (sum allfs)
        allps (reduce (fn[m [k v]] (assoc m k (float (/ v sz)))) {} allfs)]
    [ccfs&ps allfs allps sz]))


(defn cc-freqn
  "cc-freqs for freqn - use freqn over colls"
  [n colls & {par :par :or {par 1}}]
  (cc-freqs n colls freqn :par par))

(defn cc-tfreqn
  "cc-tfreqs for freqn, cc-freqs with freqn plus overall totals."
  [n colls & {par :par :or {par 1}}]
  (cc-tfreqs n colls freqn :par par))

(defn cc-freqs-probs
  "cc-freqs&probs using freqn as base frequency function"
  [n colls & {par :par :or {par 1}}]
  (cc-freqs&probs n colls #(freqs&probs %1 %2 freqn) :par par))


(defn combin-count-reduction [item-coll sym]
  (coalesce-xy-yx (freqn 1 item-coll)
                  (fn[x v] (if (not v) 0 (+ (val x) v)))))


(defn combins-freqn
  "Compute frequency map for coll where items are formed by all
   choices of coll things taken n at a time: (combins n coll).  If SYM
   is true (the default), treat vec/seq keys as symmetric.  Examples:

   (combins-freqn
     2 \"UUAAAAAACAAAAAAAAAAAACAAAACACAAAAACAAAUUCUACAAAAAAUAAAAACA\")
   {[\\C \\C] 28, [\\A \\C] 352, [\\A \\A] 946, [\\U \\C] 48,
    [\\U \\A] 264, [\\U \\U] 15}

   (combins-freqn
     2 \"UUAAAAAACAAAAAAAAAAAACAAAACACAAAAACAAAUUCUACAAAAAAUAAAAACA\" :sym false
   {[\\C \\U] 23, [\\C \\C] 28, [\\C \\A] 149, [\\A \\U] 131,
    [\\A \\C] 203, [\\A \\A] 946, [\\U \\C] 25, [\\U \\A] 133, [\\U \\U] 15}
  "
  [n coll &
   {sym :sym :or {sym true}}]
  (combin-count-reduction (combins n coll) sym))

(defn choose-k-freqn
  "Synonym for combins-freqn"
  [n coll]
  (combins-freqn n coll))

(defn cc-combins-freqn
  "combins-freqn over all collections C in colls.  Basically cc-freqs
   with freq-fn equal combins-freqn"
  [n colls & {par :par :or {par 1}}]
  (cc-freqs n colls combins-freqn :par par))

(defn combins-freqs-probs
  "freqs&probs with freq-fn equal to combins-freqn.  So, triple with
  fs and ps maps based on combinations of coll items taken n at time,
  i.e., keys are such items, with cnt element (nCk (count coll) n)."
  [n coll]
  (freqs&probs n coll combins-freqn))

(defn cc-combins-freqs-probs
  "cc-freqs&probs with freqs&probs-fn equal to combins-freqs-probs.
   So, computes freqs and probs for each C in colls based on
   combinations of coll items taken n at a time (the keys for the
   maps).  Returns [ccfsps allfs allps cnt], where ccfsps is the seq
   of triples for each C in colls, allfs is the map for all freqs,
   allps is corresponding map for all probs and cnt is total count of
   all elements.
  "
  [n colls &
   {par :par :or {par 1}}]
  (cc-freqs&probs n colls combins-freqs-probs :par par))




(defn shannon-entropy
  "Returns the Shannon entropy of a sequence: -sum(* pi (log pi)),
   where i ranges over the unique elements of S and pi is the
   probability of i in S: (freq i s)/(count s)"
  [s & {logfn :logfn :or {logfn log2}}]
  (let [fs (frequencies s)
        cnt (double (sum fs))
        H (sum (fn[[c v]]
                 (let [p (double (/ v cnt))] (double (* p (double (logfn p))))))
               fs)]
    (double (- H))))

(defn joint-entropy
  "Returns the joint entropy of S1 and S2: -sum(* pxy (log pxy)),
   where the joint probability distribution pxy is obtained from an
   exhaustive combination of elements of S1 and S2 taken 2 at a time,
   where each pair has exactly one element from S1 and exactly one
   element from S2.
  "
  [s1 s2 & {logfn :logfn sym :sym :or {logfn log2 sym true}}]
  (let [fs (combin-count-reduction (for [x s1 y s2] [x y]) sym)
        sz (sum fs)
        ps (reduce (fn[m [k v]] (assoc m k (double (/ v sz)))) {} fs)
        entropy (fn[probs] (- (sum #(* % (logfn %)) probs)))]
    (entropy (vals ps))))

(defn seq-joint-entropy
  "Returns the joint entropy of a sequence with itself: -sum(* pi (log
   pi)), where probabilities pi are of combinations of elements of S
   taken 2 at a time.
  "
  [s & {logfn :logfn :or {logfn log2}}]
  (let [[fs ps cnt] (combins-freqs-probs 2 s)
        entropy (fn[probs] (- (sum #(* % (logfn %)) probs)))]
    (entropy (vals ps))))


(defn relative-entropy
  "Take two distributions (presumably over the same space) and compute
   the expectation of their log ratio: Let px be the PMF of pdist1 and
   py be the PMF pdist2, return

   (sum (fn[px py] (* px (log2 (/ px py)))) xs ys)

   Here, pdist(1|2) are maps giving the probability distributions (and
   implicitly the pmfs), as provided by freqs-probs, probs,
   cc-freqs-probs, combins-freqs-probs, cc-combins-freqs-probs,
   et. al.  Or any map where the values are the probabilities of the
   occurance of the keys over some sample space.
  "
  [pdist1 pdist2]
  (sum (fn[pi qi]
         (let [[_ pi] pi
               [_ qi] qi]
           (* pi (log2 (/ pi qi)))))
       :|| pdist1 pdist2))

(defn KLD "Synonym for relative-entropy"[x y] (relative-entropy x y))
(defn DX||Y "Synonym for relative-entropy" [x y] (relative-entropy x y))



(defn log-odds [frq1 frq2]
  (log2 (/ frq1 frq2)))

(defn lod-score [qij pi pj]
  (log2 (/ qij (* pi pj))))

(defn raw-lod-score [qij pi pj & {scaling :scaling :or {scaling 1.0}}]
  (if (= scaling 1.0)
    (lod-score qij pi pj)
    (int (/ (lod-score qij pi pj) scaling))))




(defn joint-prob [Omega X Y & {constr :constr :or {constr true}}]
  (let [d (nCk 8 4)]
    (for [x (range 0 4) y (range 0 3)
          :let [z (- 4 x y)]
          :when (and (<= 0 z 3)
                     (= (+ x y z) 4))]
      [x y z (double (/ (* (nCk 3 x) (nCk 2 y) (nCk 3 z)) d))])))


(defn binomial-dist [N p]
  (let [q (- 1 p)]
    (for [n (range (inc N))]
      [n (* (nCk N n) (math/expt p n) (math/expt q (- N n)))])))

(defn poisson-dist [mu]
  (for [k (range 10)]
    [k (/ (* (math/expt mu k) (math/expt Math/E (- mu)))
          (n! k))]))


;;; ----------------------------------------------------------------------
;;;


;;; Fixed cost Edit distance.
;;;
;;; (levenshtein "this" "")
;;; (assert (= 0 (levenshtein "" "")))
;;; (assert (= 3 (levenshtein "foo" "foobar")))
;;; (assert (= 3 (levenshtein "kitten" "sitting")))
;;; (assert (= 3 (levenshtein "Saturday" "Sunday")))
;;; (assert (= 22 (levenshtein
;;;   "TATATTTGGAGTTATACTATGTCTCTAAGCACTGAAGCAAA"
;;;   "TATATATTTTGGAGATGCACAT"))
;;;
(defn- new-row [prev-row row-elem t]
  (reduce
   (fn [row [d-1 d e]]
     (conj row (if (= row-elem e) d-1 (inc (min (peek row) d d-1)))))
    [(inc (first prev-row))]
    (map vector prev-row (next prev-row) t)))

(defn levenshtein
  "Compute the Levenshtein (edit) distance between S and T, where S
   and T are either sequences or strings.

   Examples:  (levenshtein [1 2 3 4] [1 1 3]) ==> 2
              (levenshtein \"abcde\" \"bcdea\")   ==> 2
  "
  [s t]
  (cond
   (or (= s t "") (and (empty? s) (empty? t))) 0
   (= 0 (count s)) (count t)
   (= 0 (count t)) (count s)
   :else
   (peek (reduce
          (fn [prev-row s-elem] (new-row prev-row s-elem t))
          (range (inc (count t)))
          s))))


(defn diff-fn
  "Return the function that is 1-F applied to its args: (1-(apply f
   args)).  Intended for normalized distance metrics.

   Ex: (let [dice-diff (diff-fn dice-coeff) ...]
         (dice-diff some-set1 some-set2))
  "
  [f]
  (fn [& args] (- 1 (apply f args))))


(defn dice-coeff [s1 s2]
  (/ (* 2 (count (set/intersection s1 s2)))
     (+ (count s1) (count s2))))

(defn jaccard-index [s1 s2]
  (/ (count (set/intersection s1 s2))
     (count (set/union s1 s2))))

(defn tversky-index
  "Tversky index of two sets S1 and S2.  A generalized NON metric
   similarity 'measure'.  Generalization is through the ALPHA and BETA
   coefficients:

   TI(S1,S2) = (/ |S1^S2| (+ |S1^S2| (* ALPHA |S1-S2|) (* BETA |S2-S1|)))

   For example, with alpha=beta=1,  TI is jaccard-index
                with alpha=beta=1/2 TI is dice-coeff
   "
  [s1 s2 alpha beta]
  (let [s1&s2 (count (set/intersection s1 s2))
        s1-s2 (count (set/difference s1 s2))
        s2-s1 (count (set/difference s2 s1))]
    (/ s1&s2
       (+ s1&s2 (* alpha s1-s2) (* beta s2-s1)))))


(def
 ^{:doc
   "Named version of (diff-fn jaccard-index s1 s2).  This difference
    function is a similarity that is a proper _distance_ metric (hence
    usable in metric trees like bk-trees)."
   :arglists '([s1 s2])}
 jaccard-dist
 (diff-fn jaccard-index))


(defn freq-jaccard-index
  ""
  [s1 s2]
  (let [freq-s1 (set s1)
        freq-s2 (set s2)
        c1 (sum (set/intersection freq-s1 freq-s2))
        c2 (sum (set/union freq-s1 freq-s2))]
    (/ c1 c2)))


(defn bi-tri-grams [s]
  (let [bi-grams (set (keys (freqn 2 s)))
        tri-grams (set (keys (freqn 3 s)))]
    [(set/union bi-grams tri-grams)
     [bi-grams tri-grams]]))

(defn all-grams [s]
  (let [all-gram-sets
        (for [n (range 1 (count s))]
          (-> (freqn n s) keys set))]
    [(apply set/union all-gram-sets) all-gram-sets]))

(defn ngram-compare
  ""
  [s1 s2 & {uc? :uc? n :n scfn :scfn ngfn :ngfn
            :or {n 2 uc? false scfn dice-coeff ngfn word-letter-pairs}}]
  (let [s1 (ngfn n (if uc? (str/upper-case s1) s1))
        s2 (ngfn n (if uc? (str/upper-case s2) s2))]
    (scfn s1 s2)))

;;;(float (ngram-compare "FRANCE" "french"))
;;;(float (ngram-compare "FRANCE" "REPUBLIC OF FRANCE"))
;;;(float (ngram-compare "FRENCH REPUBLIC" "republic of france"))
;;;(float (ngram-compare
;;;        "TATATTTGGAGTTATACTATGTCTCTAAGCACTGAAGCAAA"
;;;        "TATATATTTTGGAGATGCACAT"))


(defn normed-codepoints [s]
  (vec (map #(let [nc (- % 97)]
               (cond
                (>= nc 0) nc
                (= % 32) 27
                :else 28))
            (str/codepoints s))))

(defn ngram-vec [s & {n :n :or {n 2}}]
  (let [ngrams (word-letter-pairs s n)
        ngram-points (map (fn [[x y]]
                            (int (+ (* x 27) y)))
                          (map normed-codepoints ngrams))
        v (int-array 784 0)]
    (doseq [i ngram-points]
      (aset v i 1))
    v))

