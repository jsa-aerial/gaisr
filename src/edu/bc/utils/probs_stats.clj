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

(defn keysort [map]
  (sort #(string-less? (key %1) (key %2)) map))


(defn- freqm-probm
  "Helper for some functions that can't make direct use of other
   higher level capabilities.  Mostly (only?) needed when there needs
   to be a choice between symmetric and non symmetric frequency
   distributions.
  "
  [freq-map]
  (let [fs freq-map
        sz (double (sum fs))]
    (reduce (fn[m [k v]] (assoc m k (double (/ v sz))))
            {} fs)))


(defn freqn
  "Frequencies of n-grams (aka l-mers \"features\", et. al.) in
   collection COLL treated as a sequence.  n is the \"window\" or
   resolution width.  Slide is always fixed at 1 (one) position.

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
      (let [k (apply str (take n s))]
        (if (>= (count s) n)
          (recur (rest s)
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
        sz (double (sum (vals freqs)))
        probs (reduce (fn[m [k v]]
                        (assoc m k (double (/ v sz))))
                      {} freqs)]
    [freqs probs sz]))

(defn freqs-probs
  "freqs&probs for n and coll with freq-fn equal freqn"
  [n coll]
  (freqs&probs n coll freqn))

(defn probs
  "Probabilities for items from coll taken n at a time (see freqn) or
   as determined by the frequency dist map freq-dist-map."
  ([n coll]
     (second (freqs-probs n coll)))
  ([freq-dist-map]
     (second (freqs&probs 0 freq-dist-map (fn[n coll] coll)))))


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
        allps (reduce (fn[m [k v]] (assoc m k (double (/ v sz)))) {} allfs)]
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


(defn combin-count-reduction
  "Performs (freqn 1 item-coll), giving result M.  If sym is true,
   coalesces items in M whose keys are reversible by retaining one key
   with the sum of the values of previouis two.  See coalesce-xy-yx
   for complete description.  This was originally a one off
   specifically for coalescing combination sets (see combins), hence
   the (now bogus) name.
  "
  [item-coll sym]
  (let [fm (freqn 1 item-coll)]
    (if sym
      (coalesce-xy-yx fm (fn[x v] (if (not v) 0 (+ (val x) v))))
      fm)))


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


(defn joint-prob-x
  "Xperimental.  Not sure how to make this work in general while still
   being useful.  Intent would be to supply sample space, random
   variables X & Y (as functions - including maps & vectors), that
   obey and/or define a constraint CONSTR, which may involve further
   random variables.  Generate all joint occurances and from there the
   joint probabilities.
  "
  [Omega X Y & {constr :constr :or {constr true}}]
  (let [d (nCk 8 4)]
    (for [x (range 0 4) y (range 0 3)
          :let [z (- 4 x y)]
          :when (and (<= 0 z 3)
                     (= (+ x y z) 4))]
      [x y z (double (/ (* (nCk 3 x) (nCk 2 y) (nCk 3 z)) d))])))

(defn joint-probability
  "Given a set of collections c1, c2, c3, .. cn, and combinator, a
   function of n variables which generates joint occurances from {ci},
   returns the joint probability distribution.  If sym? is true
   coalesce reversable element occurances.

   Ex: (apply joint-probability
              transpose true (take 2 (seq/rotations \"GGCGGAAGACCGCCUCGA\")))
   => {[\\U \\C] 0.11111111, [\\A \\G] 0.2777778, [\\C \\G] 0.2777778,
       [\\A \\C] 0.055555556, [\\A \\A] 0.055555556, [\\G \\G] 0.11111111,
       [\\C \\C] 0.11111111}

   Ex: (joint-probability #(combins 2 %) true \"GGCGGAAGACCGCCUCGA\")
   => {[\\C \\C] 0.09803921568627451, [\\G \\G] 0.13725490196078433,
       [\\A \\A] 0.0392156862745098, [\\A \\C] 0.1568627450980392,
       [\\C \\G] 0.27450980392156865, [\\A \\G] 0.1830065359477124,
       [\\U \\A] 0.026143790849673203, [\\G \\U] 0.0457516339869281,
       [\\U \\C] 0.0392156862745098}

   Ex: (joint-probability
         (fn[X Y]
           (for [x X]
             (if (even? x)
                 [:e (.charAt Y 0)]
                 [(second (div x 3)) (rand-nth (rest Y))])))
         false
         (take 100 (iterate inc 0)) \"abcde\")
   => {[2 \\b] 0.04, [1 \\b] 0.05, [2 \\c] 0.04, [0 \\b] 0.07, [1 \\c] 0.05,
       [2 \\d] 0.04, [0 \\c] 0.05, [1 \\d] 0.05, [2 \\e] 0.04, [0 \\d] 0.04,
       [1 \\e] 0.02, [0 \\e] 0.01, [:e \\a] 0.5}
   "
  ([combinator sym? coll]
     (freqm-probm (combin-count-reduction (combinator coll) sym?)))
  ([combinator sym? coll1 coll2]
     (freqm-probm (combin-count-reduction (combinator coll1 coll2) sym?)))
  ([combinator sym? coll1 coll2 & colls]
     (freqm-probm (combin-count-reduction
                   (apply combinator (conj colls coll1 coll2)) sym?))))

(defn JPSxy
  "Joint Probability Symmetric: joint-probability with sym? true"
  ([combinator coll]
     (joint-probability combinator true coll))
  ([combinator coll & colls]
     (apply joint-probability combinator true (cons coll colls))))

(defn JPxy
  "Joint Probability non symmetric: joint-probability with sym? false"
  ([combinator coll]
     (joint-probability combinator true coll))
  ([combinator coll & colls]
     (apply joint-probability combinator true (cons coll colls))))


(defn cond-probability
  "Given collections coll1 and coll2, and combinator, a function of 2
   variables which generates joint occurances from coll1 & coll2,
   returns the conditional probability distributions for coll1
   conditioned on elements of coll2.  If sym? is true coalesce
   reversable element occurances during joint occurance
   computation.

   Alternatively, in the distribution version, PXY is a joint
   distribution in map form, and PY is a distribution in map form.
   That is, in this variant, the probability distributions are
   explicitly provided.

   OCCURS? is a function of two variables, pkx and pkxy, where pkx is
   a key derived from coll2/PY and pkxy is a key derived from the
   combinator map of coll1 with coll2 or a key in PXY.  It returns
   whether pkx occurs in pkxy.

   cond-probability returns a map of maps, where each inner map is a
   distribution of coll1 given ei in coll2:

   {{ei {P(COLL1)|ei}} ei in COLL2.

   Use pXY|y to obtain any given distribution, where y = some ei.

   Examples:

   (cond-probability
     transpose false #(= %1 (.charAt %2 1)) \"defdeef\" \"abcaabb\")
   => {\\a {\"ea\" 0.3333333333333333, \"da\" 0.6666666666666666},
       \\b {\"fb\" 0.3333333333333333, \"eb\" 0.6666666666666666},
       \\c {\"fc\" 0.9999999999999997}}

   (cond-probability
     transpose false #(= %1 (.charAt %2 1)) \"abcaabb\" \"defdeef\")
   => {\\d {\"ad\" 1.0},
       \\e {\"ae\" 0.3333333333333333, \"be\" 0.6666666666666666},
       \\f {\"bf\" 0.5, \"cf\" 0.5}}
  "
  ([PXY PY occurs?]
     (reducem
      (fn[[pkx pvx] [pkxy pvxy]]
        #_(prn :f pkx pvx pkxy pvxy)
        (if (occurs? pkx pkxy)
          {pkx {pkxy  (/ pvxy pvx)}}
          {}))
      (fn[m subm]
        #_(prn :FR m subm)
        (merge-with into m subm))
      PY PXY))
  ([combinator sym? occurs? coll1 coll2]
     (let [PXY (joint-probability combinator sym? coll1 coll2)
           PY  (probs 1 coll2)]
       (cond-probability PXY PY occurs?))))


(defn pXY|y
  "Given a total conditional distribution XY (see cond-probability),
   return the distribution of XY with respect to y.  XY should be a
   map of maps of the form returned from cond-probability.  Typically,
   this would be used in conjunction with cond-probability, but need
   not be as long as the form of the map of maps is the same, where
   each inner map is the joint distribution conditioned on an element.

   Ex:

   (pXY|y (cond-probability transpose false #(= %1 (.charAt %2 1))
                            \"defdeef\" \"abcaabb\")
          \\a)
  => {\"ea\" 0.3333333333333333, \"da\" 0.6666666666666666}
  "
  [XY y] (XY y))


(defn binomial-pdf
  "Return the binomial probability mass/distribution _function_
   corresponding to P, the probability of a hit in some Bernoulli
   trial.  The resulting function returned, call it bnpdf, is a
   function of two parameters, the number of hits k (successes), and
   the number of trials N.  So, if the result of (binomial-pdf 1/4) is
   bound to bnpdf, then (bnpdf 3 10) would be the probability of
   getting 3 hits in 10 trials for sample probability 1/4.  A
   distribution of bnpdf can then be generated via some technique
   like (map #(bnpdf % 1000) (range 1 (inc 1000))).
  "
  [p]
  (let [q (- 1 (double p))]
    (fn[k N] (* (nCk N k) (math/expt p k) (math/expt q (- N k))))))

(defn binomial-dist
  "Binomial probability distribution. P is the probability of a hit
   and N is the number of trials.  Returns a seq of pairs [k pr],
   where k in (range 1 (inc N)) and pr is probability of getting k
   successes in N trials.
  "
  [N p]
  (if (= N 0)
    0.0
    (let [q (- 1 (double p))]
      (for [k (range 1 (inc N))]
        [k (* (nCk N k) (math/expt p k) (math/expt q (- N k)))]))))


(defn poisson-pdf
  "Return the Poisson probability mass/distribution _function_
   corresponding to mu, the mean or expected value of a random
   variable X over some time/space interval. In the Poisson case the
   mean is also the variance of X.

   The resulting function returned, call it pspdf, is a function of
   one parameter, the number of occurances K.  NOTE: pspdf is only
   defined for integer values.

   So, if the result of (poisson-pdf 1.0) is bound to pspdf,
   then (pspdf 3) would be the probability of getting 3 hits.  A
   distribution of pspdf can then be generated via some technique
   like: (map pspdf (range 1 (inc 100))).
  "
  [mu]
  (fn[k] {:pre [(integer? k)]}
    (/ (* (math/expt mu k) (math/expt Math/E (- mu)))
       (n! k))))

(defn poisson-dist
  "Poisson probability distribution. "
  [N mu]
  (for [k (range N)] ; effectively k-1
    [k (/ (* (math/expt mu k) (math/expt Math/E (- mu)))
          (n! k))]))

(defn poisson-cdf
  "Return the Poisson cumulative probability _function_ corresponding
   to mu, the mean or expected value of a random variable X over some
   time/space interval.  In the Poisson case, the mean is also the
   variance of X.

   The resulting function returned, call it pscdf, is a function of
   two parameters: lower limit l, and upper limit u, both >= 0.  l and
   u-1 are the minimum and maximum number of occurances, i.e., the
   interval is the half open [0, u) interval.

   Ex: If 2 bugs/week are introduced on average in statware, what is
   probability there will be less than 5 bugs introduced next month.
   So, mu = 2/week * 4 week = 8, l = 0 and u = 5:

   ((poisson-cdf 8) 0 5) => 0.0996 or about 10% chance.  Of course,
   those dudes aren't using Clojure! :)
  "
  [mu]
  (fn[l u] {:pre [(integer? l) (integer? u) (>= u l 0)]}
    (sum (map (poisson-pdf mu) (range l u)))))

(defn poisson-sample
  "Generate a \"random\" sample from a Poisson distribution
   characterized by mean/variance mu.  In the two parameter case,
   generate N such samples.
  "
  ([mu]
     (let [L (math/expt Math/E (- mu))]
       (loop [k 1 p 1.0]
         (if (<= p L)
           (- k 1)
           (recur (inc k) (* p (rand)))))))
  ([N mu]
     (for [k (range (inc N))]
       (poisson-sample mu))))


(defn geometric-pdf
  "Return the geometric probability mass/distribution _function_
   corresponding to p, the probability of a hit in some Bernoulli
   trial.

   The resulting function returned, call it gmpdf, is a function of
   one parameter, the number of occurances K.  NOTE: gmpdf is only
   defined for integer values.

   So, if the result of (geometric-pdf 1/4) is bound to gmpdf,
   then (gmpdf 3) would be the probability of getting a hit after 3
   tries.  A distribution of gmpdf can then be generated via some
   technique like:

   (map pspdf (range 10))
  "
  [p]
  (let [p (double p)
        q (- 1.0 p)]
    (fn[k] (* (math/expt q k) p))))

(defn geometric-dist
  "Geometric probability distribution.  P is the probability of an
   event.  Return seq of pairs [k pr], where k in (range N) and pr is
   probability of success after k failures.
  "
  [N p]
  (let [p (double p)
        q (- 1.0  p)]
    (for [k (range N)]
      [k (* (math/expt q k) p)])))




(defn pdsum
  ""
  [f pdf1 pdf2 & pdfs]
  (let [parallel (= pdf1 :||)
        pdfs (cons pdf2 pdfs)
        pdfs (if parallel pdfs (cons pdf1 pdfs))
        Omega (apply set/union (map #(set (keys %)) pdfs))]
    (apply sum (fn[k]
                  (let [pdf-vals (map #(get % k (double 0.0)) pdfs)]
                    (apply f pdf-vals)))
           pdfs)))



(defn shannon-entropy
  "Returns the Shannon entropy of a sequence: -sum(* pi (log pi)),
   where i ranges over the unique elements of S and pi is the
   probability of i in S: (freq i s)/(count s)"
  [s & {logfn :logfn :or {logfn log2}}]
  (let [fs (frequencies s)
        cnt (double (sum fs))
        H (sum (fn[[_ v]]
                 (let [p (double (/ v cnt))]
                   (double (* p (double (logfn p))))))
               fs)]
    (double (- H))))

(defn entropy
  "Entropy calculation for the probability distribution dist.
   Typically dist is a map giving the PMF of some sample space.  If it
   is a string or vector, this calls shannon-entropy on dist.
  "
  [dist & {logfn :logfn :or {logfn log2}}]
  (if (or (string? dist) (vector? dist))
    (shannon-entropy dist)
    (let [dist (if (map? dist) (vals dist) dist)]
      (- (sum (fn[v]
                (let [v (double v)]
                  (if (= 0.0 v)
                    0.0
                    (double (* v (logfn v))))))
              dist)))))


(defn joint-entropy
  "Given a set of collections c1, c2, c3, .. cn, and combinator, a
   function of n variables which generates joint occurances from {ci},
   returns the joint entropy over all the set:

   -sum(* px1..xn (log2 px1..xn))

   If sym?, treat [x y] and [y x] as equal.
  "
  ([combinator sym? coll]
     (entropy (joint-probability combinator sym? coll)))
  ([combinator sym? coll1 coll2]
     (entropy (joint-probability combinator sym? coll1 coll2)))
  ([combinator sym? coll1 coll2 & colls]
     (entropy (apply joint-probability combinator sym? coll1 coll2 colls))))

(defn HXY
  "Synonym for joint-entropy"
  [& args]
  (apply joint-entropy args))


(defn cond-entropy
  "Given the joint probability distribution PXY and the distribution
   PY, return the conditional entropy of X given Y = H(X,Y) - H(Y).

   Alternatively, given a set of collections c1, c2, c3, .. cn, and
   combinator, a function of n variables which generates joint
   occurances from {ci}, returns the multivariate conditional entropy
   induced from the joint probability distributions.  If sym? is true
   coalesce reversable element occurances in these distributions.

   For the case of i > 2, uses the recursive chain rule for
   conditional entropy (bottoming out in the two collection case):

   H(X1,..Xn-1|Xn) = (sum H(Xi|Xn,X1..Xi-1) (range 1 (inc n)))

  "
  ([PXY PY]
     (- (entropy PXY) (entropy PY)))
  ([combinator sym? coll1 coll2]
     (- (entropy (joint-probability combinator sym? coll1 coll2))
        (entropy (probs 1 coll2))))
  ([combinator sym? coll1 coll2 & colls]
     ;; Apply chain rule...
     (- (entropy (apply joint-probability combinator sym? coll1 coll2 colls))
        (apply cond-entropy combinator sym? coll2 colls))))

(defn HX|Y
  "Synonym for cond-entropy"
  [& args]
  (apply cond-entropy args))


(defn combin-joint-entropy
  "Given a set of collections c1, c2, c3, .. cn, return the joint
   entropy: - (sum (* px1..xn (log2 px1..xn))) all-pairs-over {ci}.
   Where all-pairs-over is an exhaustive combination of elements of
   {ci} taken n at a time, where each n-tuple has exactly one element
   from each ci (i.e., the cross product of {ci})
  "
  ([coll1 coll2]
     (joint-entropy
      (fn[s1 s2] (reducem (fn[x y] [[x y]]) concat s1 s2)) true coll1 coll2))
  ([coll1 coll2 & colls]
     (apply joint-entropy
            (fn[& cs]
              (apply reducem (fn[& args] [args]) concat cs))
            true coll1 coll2 colls)))


(defn seq-joint-entropy
  "Returns the joint entropy of a sequence with itself: -sum(* pi (log
   pi)), where probabilities pi are of combinations of elements of S
   taken 2 at a time.  If sym?, treat [x y] and [y x] as equal.
  "
  [s & {:keys [sym? logfn] :or {sym? false logfn log2}}]
  (joint-entropy #(combins 2 %) sym? s))


(defn relative-entropy
  "Take two distributions (that must be over the same space) and
   compute the expectation of their log ratio: Let px be the PMF of
   pdist1 and py be the PMF pdist2, return

   (sum (fn[px py] (* px (log2 (/ px py)))) xs ys)

   Here, pdist(1|2) are maps giving the probability distributions (and
   implicitly the pmfs), as provided by freqs-probs, probs,
   cc-freqs-probs, combins-freqs-probs, cc-combins-freqs-probs,
   et. al.  Or any map where the values are the probabilities of the
   occurance of the keys over some sample space.  Any summation term
   where (or (= px 0) (= py )), is taken as 0.0.

   NOTE: maps should have same keys! If this is violated it is likely
   you will get a :negRE exception or worse, bogus results.  However,
   as long as the maps reflect distributions _over the same sample
   space_, they do not need to be a complete sampling (a key/value for
   all sample space items) - missing keys will be included as 0.0
   values.

   Also known as Kullback-Leibler Divergence (KLD)

   KLD >= 0.0 in all cases.
  "
  [pdist1 pdist2]
  ;; Actually, we try to 'normalize' maps to same sample space.
  (let [Omega (set/union (set (keys pdist1)) (set (keys pdist2)))
        KLD (sum (fn[k]
                   (let [pi (get pdist1 k 0.0)
                         qi (get pdist2 k 0.0)]
                     (if (or (= pi 0.0) (= qi 0.0))
                       (double 0.0)
                       (* pi (log2 (/ pi qi))))))
                 Omega)]
    (cond
     (>= KLD 0.0) KLD
     (< (math/abs KLD) 1.0E-10) 0.0
     :else
     (raise
      :type :negRE :KLD KLD :Pdist pdist1 :Qdist pdist2))))


(defn KLD "Synonym for relative-entropy"[x y] (relative-entropy x y))
(defn DX||Y "Synonym for relative-entropy" [x y] (relative-entropy x y))


(defn mutual-information
  "Mutual information between the content of coll1 and coll2 as
   combined by combinator, a function of two variables, presumed to be
   over coll1 and coll2 returning a collection of combined elements.
   If sym? is true, treat elements with reverses as equal.

   Mutual information is the relative entropy (KLD) of the joint
   probability distribution to the product distribution.  Here the
   distributions arise out of the frequencies over coll1 and coll2
   individually and jointly over the result of combinator on coll1 and
   coll2.

   Let C be (combinator coll1 coll2).  Computes

   (sum (* pxy (log2 (/ pxy (* px py)))) xs ys) =

   (+ (shannon-entropy coll1) (shannon-entropy coll2) (- (joint-entropy C))) =

   Hx + Hy - Hxy = I(X,Y)

   (<= 0.0 I(X,Y) (min [Hx Hy]))
  "
  [combinator sym? coll1 coll2]
  (let [HX  (shannon-entropy coll1)
        HY  (shannon-entropy coll2)
        HXY (joint-entropy combinator sym? coll1 coll2)
        IXY (+ HX HY (- HXY))]
    (cond
     (>= IXY 0.0) IXY
     (< (math/abs IXY) 1.0E-10) 0.0
     :else
     (raise
      :type :negIxy :Ixy IXY :Hxy HXY :Hx HX :Hy HY))))

(defn IXY "Synonym for mutual information"
  [combinator sym? coll1 coll2]
  (mutual-information combinator sym? coll1 coll2))


(defn conditional-mutual-information
  "Conditional mutual information I(X;Y|Z).  Conditional mutual
   information is the relative entropy (KLD) of the conditional
   distribution (of two random variables on a third) with the product
   distribution of the distributed conditionals.

   Here these arise out of frequencies for the information in collz
   individually, collx & collz and colly & collz jointly as the result
   of combinator, and collx & colly & collz jointly as the result of
   combinator as well.  Hence, combinator needs to be able to handle
   the cases of 2 and 3 collections passed to it (collectively or as
   variadic). If sym? is true, treat elements with reverses as equal.

   Other interpretations/meanings: Conditional mutual information is
   the mutual information of two random variables conditioned on a
   third.  Or, put another way, it is the expected value of the mutual
   information of two RV over the values of a third.  It measures the
   amount of information shared by X & Y beyond that provided by a
   common \"mediator\".

   Let XYZ be (combinator collx colly collz)
       XZ  be (combinator collx collz)
       YZ  be (combinator colly collz)

   Computes

   sum (* pxy|z (log2 (/ pxy|z (* px|z py|z)))) xs ys zs =

    ... (some algebra and applications of Bayes) ...

   sum (* pxyz (log2 (/ (* pz pxyz) (* pxz pyz)))) xyzs xzs yzs

    ... (some more algebra and entropy identitities) ...

   (+ H(X,Z) H(Y,Z) -H(X,Y,Z) -H(Z))

   which is easier to work with...
  "
  [combinator sym? collx colly collz]
  (let [PXYZ  (vals (joint-probability combinator sym? collx colly collz))
        HXYZ  (entropy PXYZ)

        PXZ   (vals (joint-probability combinator sym? collx collz))
        HXZ   (entropy PXZ)

        PYZ   (vals (joint-probability combinator sym? colly collz))
        HYZ   (entropy PYZ)

        PZ    (probs 1 collz)
        HZ    (entropy PZ)

        IXY|Z (+ HXZ HYZ (- HXYZ) (- HZ))]
    (cond
     (>= IXY|Z 0.0) IXY|Z
     (< (math/abs IXY|Z) 1.0E-10) 0.0
     :else
     (raise
      :type :negIxy|z
      :Ixy|z IXY|Z :Hxz HXZ :Hyz HYZ :Hxyz HXYZ :Hz HZ))))

(defn IXY|Z "Synonym for conditional mutual information"
  [combinator sym? collx colly collz]
  (conditional-mutual-information combinator sym? collx colly collz))


(defn variation-information
  ""
  [combinator sym? coll1 coll2]
  (let [Hxy (joint-entropy combinator sym? coll1 coll2)
        Ixy (mutual-information combinator sym? coll1 coll2)]
    (- 1.0 (double (/ Ixy Hxy)))))


(defn total-correlation
  "One of two forms of multivariate mutual information provided here.
   The other is \"interaction information\".  Total correlation
   computes what is effectively the _total redundancy_ of the
   information in the provided content - here the information in coll1
   .. colln.  As such it can give somewhat unexpected answers in
   certain situations.

   Information content \"measure\" is based on the distributions
   arising out of the frequencies over coll1 .. colln _individually_,
   and jointly over the result of combinator applied to coll1 .. colln
   collectively.  If sym? is true, treat elements with reverses as
   equal.

   NOTE: the \"degenerate\" case of only two colls, is simply mutual
   information.

   Let C be (combinator coll1 coll2 .. colln), so xi1..xin in C is an
   element in the joint sample space, and xi in colli is an element in
   a \"marginal\" space . Computes

   sum (* px1..xn (log2 (/ px1..xn (* px1 px2 .. pxn)))) x1s x2s .. xns =

      Hx1 + Hx2 + .. + Hxn - Hx1x2..xn

   (<= 0.0
       TC(X1,..,Xn)
       (min|i (sum Hx1 .. Hxi Hxi+2 .. Hxn, i = 0..n-1, Hx0=Hxn+1=0)))

   Ex:

   (shannon-entropy \"AAAUUUGGGGCCCUUUAAA\")
   => 1.9440097497163569

   (total-correlation \"AAAUUUGGGGCCCUUUAAA\" \"AAAUUUGGGGCCCUUUAAA\")
   => 1.9440097497163569 ; not surprising

   (total-correlation
     \"AAAUUUGGGGCCCUUUAAA\" \"AAAUUUGGGGCCCUUUAAA\"
     \"AAAUUUGGGGCCCUUUAAA\" \"AAAUUUGGGGCCCUUUAAA\")
   => 5.832029249149071 ; possibly surprising if not noting tripled redundancy
  "
  ([combinator sym? coll1 coll2]
     (mutual-information combinator sym? coll1 coll2))
  ([combinator sym? coll1 coll2 & colls]
     (let [colls (cons coll1 (cons coll2 colls))
           Hxs (map shannon-entropy colls)
           HX1-Xn (apply joint-entropy combinator sym? coll1 coll2 colls)
           TC (- (sum Hxs) HX1-Xn)]
       (cond
        (>= TC 0.0) TC
        (< (math/abs TC) 1.0E-10) 0.0
        :else
        (raise
         :type :negTC
         :TC TC :Hxs Hxs :HX1-XN HX1-Xn)))))

(defn TCI
  "Synonym for total-correlation information"
  [combinator sym? coll1 coll2 & colls]
  (apply total-correlation combinator sym? coll1 coll2 colls))


(defn interaction-information
  "One of two forms of multivariate mutual information provided here.
   The other is \"total correlation\".  Interaction information
   computes the amount of information bound up in a set of variables,
   beyond that present in _any subset_ of those variables. Well, that
   seems to be the most widely held view, but there are disputes. This
   can either indicate inhibition/redundancy or facilitation/synergy
   between the variables.  It helps to look at the 3 variable case, as
   it at least lends itself to 'information' (or info venn) diagrams.

   II(X,Y,Z) = I(X,Y|Z) - I(X,Y)

   II measures the influence of Z on the information connection
   between X and Y.  Since I(X,Y|Z) < I(X,Y) is possible, II can be
   negative, as well as 0 and positive.  For example if Z completely
   determines the information content between X & Y, then I(X,Y|Z)
   would be 0 (as overlap contributed by X & Y alone is totally
   redundant), while for X & Y independent of Z, I(X,Y) could still be
   > 0.  A _negative_ interaction indicates Z inhibits (accounts for,
   causes to be irrelevant/redundant) the connection between X & Y.  A
   _positive_ interaction indicates Z promotes, is synergistic with,
   or enhances the connection.  The case of 0 is obvious...

   If we use the typical info theory venn diagrams, then II is
   represented as the area in all the circles. It is possible (and
   typically so done) to draw the Hxi (xi in {X,Y,Z}) circles so that
   all 3 have an overlap.  In this case the information between any
   two (the mutual information...) is partially accounted for by the
   third.  The extra overlap accounts for the negative value
   \"correction\" of II in this case.

   However, it is also possible to draw the circles so that any two
   overlap but all three have null intersect.  In this case, the third
   variable provides an extra connection between the other two (beyond
   their MI), something like a \"mediator\".  This extra connection
   accounts for the positive value \"facillitation\" of II in this
   case.

   It should be reasonably clear at this point that negative II is the
   much more intuitive/typical case, and that positive cases are
   surprising and typically more interesting and informative (so to
   speak...).  A positive example, would be the case of two mutually
   exclusive causes for a common effect (third variable).  In such a
   case, knowing any two completely determines the other.

   All of this can be bumped up to N variables via typical chain rule
   recurrance.

   In this implementation, the information content \"measure\" is
   based on the distributions arising out of the frequencies over
   coll1 .. colln _individually_, and jointly over the result of
   combinator applied collectively to all subsets of coll1 .. colln.
   If sym? is true, treat elements with reverses as equal.

   For collections coll1..colln and variadic combinator over any
   subset of collections, Computes:

   (sum (fn[ss] (* (expt -1 (- n |ss|))
                   (HXY combinator sym? ss)))
        (subsets all-colls))

   |ss| = cardinality/count of subset ss.  ss is a subset of
   coll1..colln.  If sym? is true, treat reversable combined items as
   equal (see HYX/joint-entropy).

   For n=3, this reduces to (- (IXY|Z combinator sym? collx colly collz)
                               (IXY combinator sym? collx colly))

   Ex:

   (interaction-information transpose false
     \"AAAGGGUUUUAAUAUUAAAAUUU\"
     \"AAAGGGUGGGAAUAUUCCCAUUU\"
     \"AAAGGGUGGGAAUAUUCCCAUUU\")
   => -1.1673250256261127

   Ex:

   (interaction-information transpose false
     \"AAAGGGUUUUAAUAUUAAAAUUU\"
     \"AAAGGGUGGGAAUAUUCCCAUUU\"
     (reverse-compliment \"AAAGGGUGGGAAUAUUCCCAUUU\"))
   => -0.282614135781623

   Ex: Shows synergy (positive value)

   (let [X1 [0 1] ; 7 independent binary vars
         X2 [0 1]
         X3 [0 1]
         X4 [0 1]
         X5 [0 1]
         X6 [0 1]
         X7 [0 1]
         X8 [0 1]
         Y1 [X1 X2 X3 X4 X5 X6 X7]
         Y2          [X4 X5 X6 X7]
         Y3             [X5 X6 X7]
         Y4                   [X7]
         bitval #(loop [bits (reverse %) n 0 v 0]
                   (if (empty? bits)
                     v
                     (recur (rest bits)
                            (inc n)
                            (+ v (* (first bits) (expt 2 n))))))
         rfn #(map bitval (apply reducem (fn[& args] [args]) concat %))
         ;; Turn to (partially redundant) value sets - each Yi is the set
         ;; of all numbers for all (count Yi) bit patterns.  So, Y1, Y2,
         ;; and Y3 share all values of 3 bits, while the group with Y4
         ;; added shares just all 1 bit patterns
         Y1 (rfn Y1)
         Y2 (rfn Y2)
         Y3 (rfn Y3)
         Y4 (rfn Y4)]
     [(interaction-information concat false Y1 Y2 Y3)
      (interaction-information concat false Y1 Y2 Y3 Y4)
      (interaction-information concat false Y1 Y2 Y4)
      (interaction-information concat false Y2 Y3 Y4)])
  => [-3.056592722891465   ; Intuitive, since 3 way sharing
       1.0803297840536592] ; Unintuitive, since 1 way sharing induces synergy!?
  "
  #_([combinator sym? collx colly collz]
       (let [Ixy|z (IXY|Z combinator sym? collx colly collz)
             Ixy   (IXY combinator sym? collx colly)]
         (- Ixy|z Ixy)))
  ([combinator sym? collx colly collz & colls]
     (let [colls (->> colls (cons collz) (cons colly) (cons collx))
           ssets (remove empty? (subsets colls))
           tcnt (count colls)]
       (- (sum (fn[ss]
                 (let [sscnt (count ss)]
                   (* (math/expt -1.0 (- tcnt sscnt))
                      (apply joint-entropy combinator sym? ss))))
               ssets)))))

(defn II
  "Synonym for interaction-information"
  [combinator sym? collx colly collz & colls]
  (apply interaction-information combinator sym? collx colly collz colls))


(defn lambda-divergence
  "Computes a symmetrized KLD variant based on a probability
   parameter, typically notated lambda, in [0..1] for each
   distribution:

   (+ (* lambda (DX||Y Pdist M)) (* (- 1 lambda) (DX||Y Qdist M)))

   Where M = (+ (* lambda Pdist) (* (- 1 lambda) Qdist))
           = (merge-with (fn[pi qi] (+ (* lambda pi) (* (- 1 lambda) qi)))
                         Pdist Qdist)

   For lambda = 1/2, this reduces to

   M = 1/2 (merge-with (fn[pi qi] (+ pi qi)) P Q)

   and (/ (+ (DX||Y Pdist M) (DX||Y Qdist M)) 2) = jensen shannon
  "
  [lambda Pdist Qdist]
  (let [Omega (set/union (set (keys Pdist)) (set (keys Qdist)))
        M (reduce (fn[M k]
                    (assoc M k (+ (* lambda (get Pdist k 0.0))
                                  (* (- 1 lambda) (get Qdist k 0.0)))))
                  {}
                  Omega)]
    (+ (* lambda (DX||Y Pdist M)) (* (- 1 lambda) (DX||Y Qdist M)))))

(defn DLX||Y
  "Synonym for lambda divergence"
  [l Pdist Qdist]
  (lambda-divergence l Pdist Qdist))


(defn jensen-shannon
  "Computes Jensen-Shannon Divergence of the two distributions Pdist
  and Qdist.  Pdist and Qdist _must_ be over the same sample space!"
  [Pdist Qdist]
  #_(let [Omega (set/union (set (keys Pdist)) (set (keys Qdist)))
          M (reduce (fn[M k]
                      (assoc M k (/ (+ (get Pdist k 0.0)
                                       (get Qdist k 0.0))
                                    2.0)))
                    {}
                    Omega)]
      (/ (+ (DX||Y Pdist M) (DX||Y Qdist M)) 2.0))
  (lambda-divergence 0.5 Pdist Qdist))




(defn log-odds [frq1 frq2]
  (log2 (/ frq1 frq2)))

(defn lod-score [qij pi pj]
  (log2 (/ qij (* pi pj))))

(defn raw-lod-score [qij pi pj & {scaling :scaling :or {scaling 1.0}}]
  (if (= scaling 1.0)
    (lod-score qij pi pj)
    (int (/ (lod-score qij pi pj) scaling))))




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

