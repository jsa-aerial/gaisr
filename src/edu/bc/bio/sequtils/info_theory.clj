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
        edu.bc.bio.seq-utils
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

(defn degap-freqs
  "Take a frequency map, and remove all elements whose keys contain
   gap characters.  Gap chars are either the defaults \\. \\- or those
   in gap-chars (a seqable collection).
  "
  ([freq-map]
     (reduce (fn[m [k v]]
               (if (gaps? k) m (assoc m k v)))
             {} freq-map))
  ([freq-map gap-chars]
     (reduce (fn[m [k v]]
               (if (gaps? k gap-chars) m (assoc m k v)))
             {} freq-map)))


(defn degap-tuples
  "Remove gap characters from a tuple of sequences.  Typically this is
   a pair of sequences (as for example arising 1from (combins 2
   some-seq-set)).  The degaping works for gaps in any (through all)
   of the elements and preserves the correct bases and their order in
   cases where gaps line up with non gaps.  Gap chars are either the
   defaults \\. \\- or those in gap-chars (a seqable collection).

   EX:

   (degap-tuples
     [\"CAAAUAAAAUAUAAUUUUAUAAUAAUAAGAAUAUAUAAAAAAUAUUAUAUAAAAGAAA\"
      \"GGGAGGGGGGGGGGG-GGGGG-GGAGGGGGGG--GGGG-GGGGGAGG-GGGG-GGGG-\"])
  => (\"CAAAUAAAAUAUAAUUUAUAUAAUAAGAAUAUAAAAAUAUUAAUAAAGAA\"
      \"GGGAGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGAGGGGGGGGGG\")
  "
  ([tuple-of-sqs]
     (transpose (filter #(not (gaps? %)) (transpose tuple-of-sqs))))
  ([tuple-of-sqs gap-chars]
     (transpose (filter #(not (gaps? % gap-chars)) (transpose tuple-of-sqs)))))




(defn- _adjust-one
  "Helper for (seqs|aln)-freqs-probs.  Basically degaps and
   recalculates freqs and probs of the passed in freqs fs and probs
   ps.
  "
  [[fs ps cnt] nogaps]
  (let [degap (if (or (string?  nogaps) (coll? nogaps))
                #(degap-freqs % nogaps)
                degap-freqs)
        reprob (fn[sz fs]
                 (reduce (fn[m [k v]] (assoc m k (double (/ v sz))))
                         {} fs))]
    (let [fs (degap fs)
          sz (sum fs)
          ps (reprob sz fs)]
      [fs ps sz])))

(defn- _str-em
  "Helper function for (seqs|aln)-freqs-probs.  Stringify the keys in
   the frequency and probability maps.
  "
  [[fs ps cnt]]
  (let [rfn #(reduce (fn[m [k v]]
                       (assoc m (apply str (ensure-vec k)) v))
                     {} %)]
    [(rfn fs) (rfn ps) cnt]))

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

   ccfs&ps is a seq of triples [fs ps cnt], for each C in seqset,
   allfs is the map of freqs over all Cs,
   allps is the map of probs over all Cs and
   tcount is the total items over coll.

   NOTE: item keys are \"stringified\", i.e., if k is a key from a map
   produced by fsps-fn, ensures that all returned maps use (apply str
   k) for all keys.
  "
  [n seqset & {:keys [fsps-fn nogaps norm par]
               :or {fsps-fn cc-freqs-probs nogaps true norm true par 1}}]
  {:pre [(gaisr-seq-set? seqset)]}
  (let [seqset (if (coll? seqset) seqset (read-seqs seqset))
        seqset (if norm (norm-elements seqset) seqset)

        [ccfsps fs ps cnt] (if (= par 1) ; Accomodate fsps-fns w/o par arg!
                             (fsps-fn n seqset)
                             (fsps-fn n seqset :par par))
        [fs ps cnt] (_str-em [fs ps cnt])
        ccfsps (pxmap _str-em par ccfsps)]
    (if (not nogaps)
      [ccfsps fs ps cnt]
      (let [ccfsps (pxmap #(_adjust-one % nogaps) par ccfsps)
            [fs ps cnt] (_adjust-one [fs ps cnt] nogaps)]
        [ccfsps fs ps cnt]))))


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
  [n seqset & {:keys [fsps-fn cols nogaps norm pgap par]
               :or {fsps-fn cc-combins-freqs-probs
                    cols false nogaps true norm true pgap 0.70 par 1}}]

  {:pre [(gaisr-seq-set? seqset)]}

  (let [aln (cond
             (and (coll? seqset) cols) (transpose seqset)
             (coll? seqset) seqset
             :else
             (read-aln-seqs seqset :cols cols))
        aln (if (= pgap 1.0) aln (filter-pgap aln pgap))

        [colsfsps fs ps cnt]
        (seqs-freqs-probs n aln :fsps-fn fsps-fn
                          :nogaps nogaps :norm norm :par par)]
    [colsfsps fs ps cnt]))


(defn bg-freqs
  "Perform a bacground frequency distribution calculation over the
   sequences in FILESPECS (a coll of legal format sequence files or a
   directory of such or if dirdir is true, a directory of directories
   of such, see read-seqs).  FTYPES gives the file types in the cases
   where filespecs is a dir or dirdir.

   By default, the distributions are performed with a sliding window
   of N length.  So, on DNA/RNA sequences 1 gives base probabilities,
   2 gives a dinucleotide distribution, etc.  To change this supply a
   different freq&prob calculation function for fsps-fn.  For more
   information see seqs-freqs-probs description.

   If cols is true, computation is over the columns of the sequences.
   if sym? is true, treat reversable keys as equal.
   If nogaps is true, removes default gap characters from calculation.
   If nogaps is a coll (for example, (keys +NONSTD-RNA+)), removes all
   those characters from calculation.  If norm is true, normalize
   characters to upper case.
  "
  [n filespecs & {:keys [fsps-fn ftypes dirdir cols sym? nogaps norm par]
                  :or {fsps-fn cc-freqs-probs ftypes [".sto"] dirdir false
                       cols false sym? false mnogaps false norm true par 1}}]
  (let [isdir (and (not (coll? filespecs))
                   (fs/directory? filespecs))
        filespecs (cond
                   (and isdir (not dirdir))
                   (apply concat
                          (map #(fs/directory-files filespecs %) ftypes))
                   dirdir filespecs ; keep as original dirdir filespec
                   :else filespecs)
        merger (fn[M m] (merge-with + M m))
        f (fn[sqs]
            (second (seqs-freqs-probs
                     n sqs :fsps-fn fsps-fn
                     :nogaps nogaps :norm norm :par par)))
        fs (cond
            (and isdir dirdir)
            (reduce (fn[M m] (merge-with + M m))
                    (map #(bg-freqs n % :fsps-fn fsps-fn
                                    :ftypes ftypes :cols cols
                                    :nogaps nogaps :norm norm :par par)
                         (fs/directory-files filespecs "")))

            (seq filespecs)
            (if cols
              (reduce-aln-seqs f merger cols filespecs)
              (reduce-seqs f merger filespecs))

            :else nil)]
    (if sym?
      (coalesce-xy-yx fs (fn [x v] (if (not v) 0 (+ (val x) v))))
      fs)))

(defn bg-freqs-probs
  "Like bg-freqs, but with the additional final computation of
   probability distribution for the frequency distribution
  "
  [n filespecs & {:keys [fsps-fn ftypes dirdir cols sym? nogaps norm par]
                  :or {fsps-fn cc-freqs-probs ftypes [".sto"] dirdir false
                       cols false sym? false nogaps false norm true par 1}}]
  (let [fqs (bg-freqs n filespecs
                      :fsps-fn fsps-fn :ftypes ftypes :dirdir dirdir
                      :cols cols :sym? sym? :nogaps nogaps :norm norm :par par)]
    [fqs (probs fqs)]))




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
  "Helper function.  Filters and transforms SQS according to switches
   COLS, NORM, and PGAP.  If cols transpose the matrix represented by
   SQS.  If norm, normalize elements of sqs by means of norm-elements.
   In all cases filter out sqs with gap content greater than pgap (a
   percentage value that defaults to 0.25)
  "
  [sqs cols norm pgap]
  (cond
   (and cols norm) (-> sqs transpose (filter-pgap pgap) norm-elements)
   cols (-> sqs transpose (filter-pgap pgap) norm-elements)
   norm (-> sqs (filter-pgap pgap) norm-elements)
   (= pgap 1.0) sqs
   :else (filter-pgap sqs pgap)))

(defn- account-for-symmetry
  "Helper function for aln-conditional-mutual-information.  Transform
   content (triples of xy base keys and associated \"residual\" column
   bases so that all xy yx keys are treated as multiple counts of xy.
   That is, treat keys as symmetrically equal.  This occurs _before_
   freq counts, so achieve this effect by literally duplicating
   elements in the resulting collection.
  "
  [triples]
  (let [x (coalesce-xy-yx triples
                          ;; just add copies of existing one
                          (fn [x v] (if (not v) () (cons (second x) v))))]
    ;; Now flatten the result - note not amenable to std flatten fn!
    (reduce (fn[c [x vs]]
              (reduce (fn[c v]
                        (conj c [x v]))
                      c vs))
            [] x)))

(defn- seq-pairs&indices
  "Helper function for *mutual-information computations.  Removes all
   seqs with more than pgap percentage of gaps and from the resulting
   set, creates the sets of seq pairs and corresponding indices.
  "
  [aln pgap par]
  (let [cnt (count aln)
        seqs&indices (partition-all 2 (interleave aln (range cnt)))
        seqs&indices (pxmap (fn[[s i]]
                              (when (seq (filter-pgap [s] pgap)) [s i]))
                            par seqs&indices)
        seqs&indices (filter #(not (empty? %)) seqs&indices) ; remove nulls
        seq-pairs&indices (map (fn[x] (transpose x)) (combins 2 seqs&indices))]
    [(map first seq-pairs&indices) (map second seq-pairs&indices)]))


(defn aln-conditional-mutual-information
  "Mutual information of all 2 column pairs in an alignment
   conditioned by the residual - unordered - bases of the remaining
   columns.  Let colpairs be (combins 2 (transpose aln)).  For any
   pair of columns [X Y] in colpairs, let Z be colpairs - {X Y}.
   Compute I(X;Y|Z), the mutual information for X&Y given
  "
  [seqset & {par :par nogaps :nogaps pgap :pgap
             cols :cols norm :norm sym? :sym?
             :or {par 1 nogaps true pgap 0.25 cols true norm true sym? true}}]
  {:pre [(gaisr-seq-set? seqset)]}
  (let [aln (if (coll? seqset) seqset (read-aln-seqs seqset))
        aln (adjust-seqs-info aln cols norm 1.0)
        [seq-pairs indices] (seq-pairs&indices aln pgap par)
        aln-map (reduce (fn[m [s1 s2]] (assoc m s1 true s2 true)) {} seq-pairs)

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
        triples (if sym? (pxmap account-for-symmetry par triples) triples)

        xyz-joint-probs (map #(probs 1 %) triples)
        HXYZ (map entropy xyz-joint-probs)

        z-probs (map (fn[tpl] (probs 1 (map second tpl))) triples)
        HZ (map entropy z-probs)

        xz-joint-probs (pxmap (fn[tpl]
                                (probs 1 (map (fn[[xy z]]
                                                [(first xy) z])
                                              tpl)))
                              par
                              triples)
        HXZ (map entropy xz-joint-probs)
        yz-joint-probs (pxmap (fn[tpl]
                                (probs 1 (map (fn[[xy z]]
                                                [(second xy) z])
                                              tpl)))
                              par
                              triples)
        HYZ (map entropy yz-joint-probs)

        Ixy|z (pxmap (fn[Hxz Hyz Hxyz Hz]
                       (let [Ixy|z (+ Hxz Hyz (- Hxyz) (- Hz))]
                         (cond
                          (>= Ixy|z 0.0) Ixy|z
                          (< (abs Ixy|z) 1.0E-10) 0.0
                          :else
                          (raise
                           :type :negIxy|z
                           :Ixy|z Ixy|z :Hxz Hxz :Hyz Hyz :Hxyz Hxyz :Hz Hz))))
                     par
                     HXZ HYZ HXYZ HZ)]
    [Ixy|z seq-pairs indices]))


(defn aln-mutual-information
  ""
  [seqset & {:keys [par nogaps pgap cols norm sym?]
             :or {par 1 nogaps true pgap 0.25 cols false norm true sym? true}}]
  {:pre [(gaisr-seq-set? seqset)]}
  (let [probs (fn[fs]
                (let [sz (sum fs)]
                  (reduce (fn[m [k v]] (assoc m k (double (/ v sz))))
                          {} fs)))
        entropy (fn[probs] (- (sum #(* % (ln %)) (vals probs))))
        aln (if (coll? seqset) seqset (read-aln-seqs seqset))
        aln (adjust-seqs-info aln cols norm 1.0)

        [seq-pairs indices] (seq-pairs&indices aln pgap par)

        pair-freqs (pxmap #(combin-count-reduction % sym?) par
                          (map transpose seq-pairs))
        pair-freqs (if nogaps (map degap-freqs pair-freqs) pair-freqs)

        joint-probs (map probs pair-freqs)
        joint-entropy (map entropy joint-probs)
        mutual-info (pxmap (fn[Hxy sp]
                             (let [[x y] (if nogaps (degap-tuples sp) sp)
                                   Hx (shannon-entropy x :logfn ln)
                                   Hy (shannon-entropy y :logfn ln)
                                   Ixy (+ Hx Hy (- Hxy))]
                               (cond
                                (>= Ixy 0.0) Ixy
                                (< (abs Ixy) 1.0E-10) 0.0
                                :else
                                (raise
                                 :type :negIxy :Ixy Ixy :Hxy Hxy :pair sp))))
                           par
                           joint-entropy seq-pairs)]
    [mutual-info joint-entropy seq-pairs pair-freqs indices]))


(defn seq-pairs-bpfreqs
  [seq-pairs & {:keys [nogaps sym? par] :or {nogaps true sym? true par 4}}]
  (pxmap (fn[p]
           (let [fm (combin-count-reduction (transpose p) sym?)]
             (if nogaps (degap-freqs fm) fm)))
         par seq-pairs))

(defn bp-stats [bp-freq-map]
  (let [sz (sum bp-freq-map)
        au (get bp-freq-map "AU" 0)
        ua (get bp-freq-map "UA" 0)
        AU (+ au ua)
        gc (get bp-freq-map "GC" 0)
        cg (get bp-freq-map "CG" 0)
        GC (+ gc cg)
        bps (+ AU GC)]
    [{:BP bps :AU AU :GC GC :au au :ua ua :gc gc :cg cg}
     {:BP (double (/ bps sz))
      :AU (double (/ AU sz))
      :GC (double (/ GC sz))}
     sz]))


(defn seq-vi
  [seqx seqy & {nogaps :nogaps :or {nogaps true}}]
  (let [seqs [seqx seqy]
        [[Ixy] [Hxy]] (aln-mutual-information seqs :pgap 1.0 :nogaps nogaps)]
    (cond
     (and (= Hxy 0.0) (not= seqx seqy)) 1
     (= Hxy 0.0) 0
     :else
     (- 1.0 (/ Ixy Hxy)))))
