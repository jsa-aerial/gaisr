;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                              U T I L S                                   ;;
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
;; Author: Jon Anthony                                                      ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns edu.bc.utils

  "General utility functions and macros.  Basically these resources
   are fairly general and intended to be usable on most any part of
   most any project"

  (:require [clojure.core.reducers :as r]
            [clojure.contrib.math :as math]
            [clojure.contrib.combinatorics :as comb]
            [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.contrib.io :as io]
            [clojure.contrib.properties :as prop]
            [edu.bc.fs :as fs])

  (:use [slingshot.slingshot
         :only [throw+ try+ get-throw-context]]
        [clojure.stacktrace
         :only [print-stack-trace]]
        [clj-shell.shell
         :only (sh)]
        [clojure.pprint
         :only [cl-format]])

  (:import (java.util Date Calendar Locale)
           java.lang.Thread
           (java.text SimpleDateFormat)))



;;; (compile 'edu.bc.utils)


;;; -----------------------------------------------------------------
;;; Miscellaneous changes/additions/"fixes" to various operations that
;;; are either already in Clojure or in bits of contribs or probably
;;; _will_ be in future (at which point these can be retired...)


(defmacro defparameter
  "def with auto declaration of symbol SYM to be dynamic"
  ([sym val]
     `(def ~(with-meta sym (assoc (or (meta sym) {}) :dynamic true)) ~val))
  ([mmap sym val]
     (let [sym-meta (or (meta sym) {})
           mmap (if (string? mmap)
                  (assoc sym-meta :doc mmap)
                  (into sym-meta mmap))]
       `(def ~(with-meta sym (assoc mmap :dynamic true)) ~val))))

(defn add-meta
  [x m]
  (with-meta x (merge (or (meta x) {}) m)))




(defn timefn
  "Time the application of F (a function) to ARGS (any set of args as
   expected by f.  Returns a two element vector [ret time] where,

   ret is the return value of f

   time is the time f took to compute ret in milliseconds
  "
  [f & args]
  (let [start-time (. java.lang.System (nanoTime))
        ret (apply f args)]
    [ret (/ (double (- (. java.lang.System (nanoTime)) start-time))
            1000000.0)]))

(defparameter {:private true} *uid* (atom (.getTime (Date.))))


(defn gen-uid
  "Generates a unique integer ID based on universal time."
  []
  (swap! *uid* inc))

(defn gen-kwuid
  "Generates a unique keyword id whose name is the str of gen-uid"
  []
  (keyword (str (gen-uid))))


(defn sleep [msecs]
  (Thread/sleep msecs))

(defn str-date
  ([] (str-date (Date.) "yyyy-MM-dd HH:mm:ss"))
  ([fm] (str-date (Date.) fm))
  ([d fm] (.format (SimpleDateFormat. fm) d)))


(defn cpu-use
  "Obtain and return cpu utilization.  Requires the 'top' command is
   available, and takes the second of two samplings (first samples
   from top are typically inaccurate).  INFO indicates the type of
   information returned:

   :idle, just the avg % idle, this is the default
   :used, just the avg % used
   :both, both idle and used

   In  all cases  uses  aggregate cpu  utilization  over SMP  systems.
   Values are returned as Doubles,  or for :both, a two element vector
   of Doubles [idle used].
  "
  [& {:keys [info] :or {info :idle}}]
  (let [idle (->> (runx "top" "-n" "2" "-b" "-d" "0.01" "-p" "1")
                  (str/split #"\n") (filter #(re-find #"^Cpu" %))
                  last (str/split #",\s+") (filter #(re-find #"%id" %))
                  first (str/split #"%") first Double.)
        use (- 100.0 idle)]
    (case info
          :idle idle
          :used use
          :both [idle use]
          idle)))


;;; Extra predicates...
(defn map-entry? [x]
  (instance? clojure.lang.MapEntry x))


(defn third [coll]
  (when (> (count coll) 2)
    (nth coll 2)))


(defn drop-until
  "Complement of drop-while"
  [pred coll] (drop-while (complement pred) coll))

(defn take-until
  "Complement of take-while"
  [pred coll] (take-while (complement pred) coll))


(defn take-until-nochange
  "Eagerly takes from SQ until consecutive elements are the same.  So,
   take until and up to element Ei, where Ei=Ei+1.  Equality of
   elements is determined by elt=, which defaults to =.
  "
  [sq & {:keys [elt=] :or {elt= =}}]
  (loop [rsq (rest sq)
         prev (first sq)
         res [prev]]
    (if (elt= prev (first rsq))
      res
      (recur (rest rsq)
             (first rsq)
             (conj res (first rsq))))))


(defn ensure-vec
  "Return a vector representation for x.  If x is a vector just return
   it, if it is a seqable return (vec x), if it is an \"atom\" return
   [x]."
  [x]
  (cond
   (vector? x) x
   (seq? x) (vec x)
   (or (map? x) (set? x)) (vec x)
   true [x]))


(defn in
  "Return whether e is an element of coll."
  [e coll]
  (if (map? coll)
    (find coll e)
    (some #(= % e) coll)))


(defn pos
  "Returns a lazy seq of positions of X within COLL taken as a sequence"
  [x coll]
  (keep-indexed #(when (= x %2) %1) coll))

(defn pos-any
  "Returns a lazy seq of positions of any element of TEST-COLL within
   COLL taken as a sequence"
  [test-coll coll]
  (keep-indexed #(when (in %2 test-coll) %1) coll))


(defn merge-with*
  "Merge-with needs to call user supplied F with the KEY as well!!!
  Returns a map that consists of the rest of the maps conj-ed onto the
  first.  If a key occurs in more than one map, the mapping(s) from
  the latter (left-to-right) will be combined with the mapping in the
  result by calling (f key val-in-result val-in-latter)."
  {:added "1.jsa"} [f & maps]
  (when (some identity maps)
    (let [merge-entry (fn [m e]
                        (let [k (key e) v (val e)]
                          (if (contains? m k)
                            (assoc m k (f k (get m k) v))
                            (assoc m k v))))
          merge2 (fn [m1 m2]
                   (reduce merge-entry (or m1 {}) (seq m2)))]
      (reduce merge2 maps))))


(defn transpose
  "Matrix transposition.  Well, effectively.  Can be used in some
   other contexts, but does the same computation.  Takes colls a
   collection of colletions, treats this as a matrix of (count colls)
   rows, each row being a string or seqable data structure: M[rij],
   where rij is the jth element of the ith row.  Returns M' = M[cji],
   where cji is the ith element of the jth column of M.
  "
  ([colls]
     {:pre [(or (string? colls)
                (and (coll? colls)
                     (every? #(or (coll? %) (string? %)) colls)))]}
     (let [colls (if (string? colls) [colls] colls)]
       (if (empty? colls)
         colls
         (let [tm (apply map vector colls)]
           (if (every? string? colls)
             (map #(apply str %) tm)
             tm)))))

  ([coll1 coll2 & colls]
     (transpose (cons coll1 (cons coll2 colls)))))


(defn- reduce-in-parallel
  "Helper function for reducem.  When reducem determines that the
   function application should proceed in parallel over the sequences,
   it defers to this function to compute the resuls"
  [f fr coll1 & colls]
  (reduce
   (fn[x y]
     (if (= x :ignore) y (fr x y)))
   :ignore (apply map f (cons coll1 colls))))

(defn reducem
  "Multiple collection reduction.  FR is a function of two arguments -
   the first is the result of the previous fold, and the second is the
   current result of F applied to the arguments from the supplied
   collections (treated as seqs).  Note, for the first application,
   the result of F is returned without invoking FR.

   By default, reduction proceeds on the results of F applied over the
   _cross product_ of the collections.  If reduction should proceed
   over collections in parallel, the first \"coll\" given should be
   the special keyword :||.  If given, this causes F to be applied to
   the elements of colls as stepped in parallel: f(coll1[i] coll2[i]
   .. colln[i]), i in [0 .. (count smallest-given-coll)].
  "
  ([f fr coll]
     (reduce
      (fn[x y]
        (if (= x :ignore)
          (f y)
          (fr x (f y))))
      :ignore coll))

  ([f fr coll1 & colls]
     (let [colls (cons coll1 colls)]
       (if (= coll1 :||)
         (apply reduce-in-parallel f fr (rest colls))
         ;;Else: We reduce X-product reductions by currying outer args
         ;;into f.
         (reduce
          (fn[r xr]
            (if (= :ignore r)
              (apply reducem (fn[& args] (apply f xr args)) fr
                     (rest colls))
              (fr r (apply reducem (fn[& args] (apply f xr args)) fr
                           (rest colls)))))
          :ignore (first colls))))))


(defn xfold
  "Fold function f over a collection or set of collections (c1, ...,
   cn) producing a collection (concrete type of vector).  Arity of f
   must be equal to the number of collections being folded with
   parameters matching the order of the given collections.  Folding
   here uses the reducer lib fold at base, and thus uses work stealing
   deque f/j to mostly relieve the partition problem.  In signatures
   with N given, N provides the packet granularity, or if (< N 1),
   granularity is determined automatically (see below) as in the
   base (non N case) signature.

   While xfold is based on r/fold, it abstracts over the combiner,
   reducer, work packet granularity, and transforming multiple
   collections for processing by f by chunking the _transpose_ of the
   collection of collections.

   As indicated above, xfold's fold combiner is monoidal on vectors:
   it constructs a new vector from an l and r vector, and has identity
   [] (empty vector).  In line with this, xfold's reducer builds up
   new vectors from elements by conjing (f a1, ... an) onto a prior
   result or the combiner identity, [], as initial result.

   Packet granularity is determined automatically (base case or N=0)
   or supplied with N > 1 in signatures with N.  Automatic
   determination tries to balance significant work chunks (keep thread
   overhead low), with chunks that are small enough to have multiple
   instances per worker queue (supporting stealing by those workers
   that finish ahead of others).
  "
  ([f coll]
     (let [cores (.. Runtime getRuntime availableProcessors)
           workers (int (math/floor (* 3/4 cores)))
           base (* 8 cores)
           n (max 2 (int (math/floor (/ (count coll) (* 2 workers)))))
           n (min base n)]
       #_(println :>>N n)
       (xfold f n coll)))
  ([f n coll]
     (if (< n 1)
       (xfold f coll)
       (r/fold n
               (fn
                 ([] [])
                 ([l r] (apply conj l r)))
               (fn[v x] (conj v (f x)))
               (vec coll))))
  ([f n coll & colls]
     (xfold (fn[v] (apply f v)) n (apply transpose coll colls))))

(defn foldx [& args] (apply xfold args))


(defn pxmap
  "Constrained pmap.  Constrain pmap to at most par threads.
   Generally, to ensure non degrading behavior, par should be
   <= (.. Runtime getRuntime availableProcessors).  It can be more,
   but if par >> availableProcessors, thrashing (excessive context
   switching) can become an issue.  Nevertheless, there are cases
   where having par be larger can reduce the ill effects of the
   partition problem.  NOTE: no effort is made to provide the true (or
   even a \"good\") solution to the partitioning of f over coll(s).

   Effectively, (pmap f (partition-all (/ (count coll) par) coll).
   Implicit doall on results to force execution.  For multiple
   collection variants, chunks the _transpose_ of the collection of
   collections.
  "
  ([f par coll]
     (if (= par 1)
       (map f coll)
       (apply concat
              (doall (pmap (fn[subset] (doall (map f subset)))
                           (partition-all (/ (count coll) par) coll))))))
  ([f par coll1 coll2]
     (if (= par 1)
       (map f coll1 coll2)
       (pxmap (fn[[x y]] (f x y)) par (transpose coll1 coll2))))
  ([f par coll1 coll2 & colls]
     (if (= par 1)
       (apply map f coll1 coll2 colls)
       (pxmap (fn[v] (apply f v)) par (apply transpose coll1 coll2 colls)))))


;;; "Had to" bring this in from COMB name space as it was private, but
;;; needed to apply it here for index return version of combins
(defn index-combinations
  "Generates all n-tuples, as indices, over set cardinality cnt"
  [n cnt]
  (lazy-seq
   (let [c (vec (cons nil (for [j (range 1 (inc n))] (+ j cnt (- (inc n)))))),
         iter-comb
         (fn iter-comb [c j]
           (if (> j n) nil
               (let [c (assoc c j (dec (c j)))]
                 (if (< (c j) j) [c (inc j)]
                     (loop [c c, j j]
                       (if (= j 1) [c j]
                           (recur (assoc c (dec j) (dec (c j))) (dec j)))))))),
         step
         (fn step [c j]
           (cons (rseq (subvec c 1 (inc n)))
                 (lazy-seq (let [next-step (iter-comb c j)]
                             (when next-step
                               (step (next-step 0) (next-step 1)))))))]
     (step c 1))))

;;; "Had to" change combins to directly use index-combinations for two
;;; return version case.
;;;
;;; Originally just (lazy-seq (map vec (comb/combinations coll k)))
(defn combins
  "Return the set of all K element _combinations_ (not permutations)
   formed from coll.  Coll does not need to be a set.  In particular,
   repeated elements are legal.

   Examples:
   (combins 2 \"abcdef\")
   => ([\\a \\b] [\\a \\c] [\\a \\d] [\\a \\e] [\\a \\f] [\\b \\c]
       [\\b \\d] [\\b \\e] [\\b \\f] [\\c \\d] [\\c \\e] [\\c \\f]
       [\\d \\e] [\\d \\f] [\\e \\f])

   (map #(apply str %) (combins 2 \"AAGGGCGUGG\"))
   => (\"AA\" \"AG\" \"AG\" \"AC\" \"AG\" \"AU\" \"AG\" \"AG\" \"AC\"
       \"AG\" \"AU\" \"GG\" \"GC\" \"GG\" \"GU\" \"GC\" \"GG\" \"GU\"
       \"CG\" \"CU\" \"GU\")
  "
  [k coll & [indices]]
  (let [v (vec (reverse coll))]
    (if (zero? k) (list [])
        (let [cnt (count coll)]
          (cond (> k cnt) nil
                (= k cnt) (list (vec coll))
                ;; Just the sets
                (not indices)
                (map #(vec (map v %)) (index-combinations k cnt))

                :else ; Sets and their indices
                [(map #(vec (map v %)) (index-combinations k cnt))
                 (combins k (range cnt))]
                )))))

(defn choose-k
  "Synonym for combins"
  [k coll]
  (combins k coll))


(defn subsets
  "All the subsets of elements of coll"
  [coll]
  (mapcat #(combins % coll) (range (inc (count coll)))))

(defn random-subset
  "Create a \"random\" N element subset of the collection s treated as a set,
   i.e., s with no duplicate elements.  If n <= 0, return #{}, the
   empty set.  If (count (set s)) <= n, return (set s).  Otherwise,
   pick N random elements from (set s) to form subset.
  "
  [s n]
  (let [s (vec (set s))]
    (cond
     (<= n 0) #{}
     (<= (count s) n) (set s)
     :else
     (loop [ss #{(rand-nth s)}]
       (if (= (count ss) n)
         ss
         (recur (conj ss (rand-nth s))))))))


(defn coalesce-xy-yx
  "Coaleseces elements of item-coll, which are or have common \"keys\",
   according to the function f.  Two keys k1 and k2 are considered
   common if (or (= k1 k2) (= (reverse k1) k2)) for reversible keys or
   simply (= k1 k2) for non reversible keys.  Reversible keys are
   vectors, seqs, or string types.

   F is a function of two parameters [x v], where x is an element of
   item-coll, and v is the current value associated with the key of x
   or nil if no association yet exists.  F is expected to return the
   current association for key of x based on x and v.  If x is a
   mapentry, (key x) is used to determine the association.  If x is a
   list or vector (first x) is used to determine the association.

   Ex:
   (freqn 1 (map #(apply str %) (combins 2 \"auauuagcgccg\")))
   => {\"aa\" 3, \"cc\" 3, \"gg\" 3, \"uu\" 3, \"ac\" 9, \"cg\" 4,
       \"ag\" 9, \"ua\" 4, \"uc\" 9, \"ug\" 9, \"au\" 5, \"gc\" 5}

   (coalesce-xy-yx *1 (fn[x v] (if (not v) 0 (+ (val x) v))))
   => {\"aa\" 3, \"cc\" 3, \"gg\" 3, \"uu\" 3, \"ac\" 9, \"cg\" 9,
       \"ag\" 9, \"ua\" 9, \"uc\" 9, \"ug\" 9}
  "
  [item-coll f]
  (let [rev (fn[x] (if (string? x) (str/reverse x) (reverse x)))
        res (reduce (fn[m x]
                      (let [k (cond
                               (map-entry? x) (key x)
                               (coll? x) (first x)
                               :else x)
                            keycoll? (or (vector? x) (string? x) (seq? x))
                            [k v] (if (not keycoll?)
                                    [k (get m k (f x nil))]
                                    (if-let [v (get m k)]
                                      [k v]
                                      (let [rk (rev k)]
                                        (if-let [v (get m rk)]
                                          [rk v]
                                          [k (f x nil)]))))]
                        (assoc m k (f x v))))
                    {} item-coll)]
    (cond
     (map? item-coll) res
     (vector? item-coll) (vec res)
     :else (seq res))))




;;; -----------------------------------------------------------------
;;; Various ancillary math/arithmetic stuff.

(defn div
  "Integer division.  Return [q r], such that floor(n / d) * q + r = q"
  [n d]
  (let [q (math/floor (/ n d))]
    [q (rem n d)]))


;;; Forward decl for sum-in-parallel
(def sum)

(defn- sum-in-parallel
  "Helper function for sum.  When sum determines that the function
   application should proceed in parallel over the sequences, it
   defers to this function to compute the results."
  ([f coll1 coll2]
     (sum (map f coll1 coll2)))
  ([f coll1 coll2 & colls]
     (let [colls (cons coll1 (cons coll2 colls))]
       (sum (apply map f colls)))))

(defn sum
  "Return the sum of the numbers in COLL.  If COLL is a map, return
   the sum of the (presumed all numbers) in (vals coll).  For function
   F versions, return the sum x in COLL (f x) or sum x in COLL1, y in
   COLL2 (f x y) or sum x1 in C1, x2 in C2, ... xn in Cn (f x1 ... xn).

   By default, summation proceeds on the results of F applied over the
   _cross product_ of the collections.  If summation should proceed
   over the collections in parallel, the first \"coll\" given should
   be the special keyword :||.  If given this causes F to be applied
   to the elements of colls as stepped in parallel: f(coll1[i]
   coll2[i] .. colln[i]), i in [0 .. (count smallest-given-coll)].

   Examples:

   (sum + [1 2 3] [1 2 3])
   => 36 ; sum of all [1 2 3] X [1 2 3] pairwise sums
   (sum + :|| [1 2 3] [1 2 3])
   => 12 ; sum of [(+ 1 1) (+ 2 2) (+ 3 3)]

   (sum (fn[x y] (* x (log2 y))) :|| [1 2 3] [1 2 3])
   => 6.754887502163469
  "
  ;; NOTE: (probably??) Should reimplement with reducem.  See for
  ;; example, prod below...
  ([coll]
     (let [vs (if (map? coll) (vals coll) coll)]
       (apply +' vs)))
  ([f coll]
     (reduce (fn [x i]
               (+' x (f i)))
             0 coll))
  ([f coll1 coll2]
     (reduce
      (fn[r c1i]
        (+' r (sum #(f c1i %1) coll2)))
      0 coll1))
  ([f coll1 coll2 & colls]
     (let [colls (cons coll1 (cons coll2 colls))]
       (if (= coll1 :||)
         (apply sum-in-parallel f (rest colls))

         ;; Else: We reduce X-product reductions by currying outer
         ;; args into f
         (reduce
          (fn[r cxi]
            (+' r (apply sum (fn[& args] (apply f cxi args)) (rest colls))))
          0 (first colls))))))


(defn cumsum
  "Cumulative sum of values in COLL.  If COLL is a map, return the
   cumulative sum of the (presumed all numbers) in (vals coll).  For
   function F versions, return the cumulative sum of x in COLL (f x)
   or the cumulative sum of x1 in COLL1, x2 in COLL2, ... xn in
   COLLn (f x1 ... xn).

   The cumulative sum is the seq of partial sums across COLL treated
   as a vector for each (i (range (count COLL))), effectively:

   [(coll 0) (sum (subvec coll 0 2)) .. (sum (subvec coll 0 (count coll)))]

   Except actual computational time is linear.
  "
  ([coll]
     (let [cv (vec (if (map? coll) (vals coll) coll))
           c0 (cv 0)]
       (first
        (reduce (fn[[v s] x]
                  (let [s (+ s x)]
                    [(conj v s) s]))
                [[c0] c0] (rest cv)))))
  ([f coll]
     (cumsum (map f coll)))
  ([f coll1 coll2 & colls]
     (cumsum (apply map f coll1 coll2 colls))))


(defn prod
  "Return the product of the numbers in COLL.  If COLL is a map,
   return the product of the (presumed all numbers) in (vals coll).
   For function F versions, return the product of x in COLL (f x) or
   prod x in COLL1, y in COLL2 (f x y) or prod x1 in C1, x2 in C2,
   ... xn in Cn (f x1 ... xn).

   By default, multiplication proceeds on the results of F applied
   over the _cross product_ of the collections.  If multiplication
   should proceed over the collections in parallel, the first \"coll\"
   given should be the special keyword :||.  If given, this causes F
   to be applied to the elements of colls as stepped in parallel:
   f(coll1[i] coll2[i] .. colln[i]), i in [0 .. (count smallest-coll)].

   Examples:

   (prod + [1 2 3] [1 2 3])
   => 172800 ; product of all [1 2 3] X [1 2 3] pairwise sums
   (prod + :|| [1 2 3] [1 2 3])
   => 48 ; product of [(+ 1 1) (+ 2 2) (+ 3 3)]
  "
  ([coll]
     (let [vs (if (map? coll) (vals coll) coll)]
       (apply *' vs)))
  ([f coll]
     (reducem f *' coll))
  ([f coll1 coll2 & colls]
     (let [colls (cons coll2 colls)]
       (apply reducem f *' coll1 colls))))


(defn sqr
  "Square x, i.e., returns (* x x)"
  [x]
  (*' x x))


(defn logb
  "Return the log to the base b _function_ of x"
  [b]
  (let [lnb (Math/log b)]
    (fn[x] (/ (Math/log x) lnb))))

(defn log
  "Return the natural log of x"
  [x]
  (Math/log x))

(defn ln
  "Return the natural log of x"
  [x]
  (Math/log x))

(def
 ^{:doc
   "Named version of (logb 2).  Important enough to have a named top
    level function"
   :arglists '([x])}
 log2 (logb 2))

(def
 ^{:doc
   "Named version of (logb 10).  Important enough to have a named top
    level function"
   :arglists '([x])}
 log10 (logb 10))


(defn n!
  "For positive integer N, compute N factorial."
  [n]
  {:pre [(integer? n) (> n -1)]}
  (reduce *' (range n 1 -1)))

(defn n-k!
  "For positive integers N and K, compute (n-k)!"
  [n k]
  {:pre [(integer? n) (integer? k) (> n -1) (<= 0 k n)]}
  (if (= n k)
    1
    (n! (-' n k))))

(defn nCk
  "For positive integers N and K, compute N choose K (binomial
   coefficient): n!/((n-k)!k!)"
  [n k]
  {:pre [(integer? n) (integer? k) (> n -1) (<= 0 k)]}
  (if (< n k)
    0
    (/ (reduce *' (range n (-' n k) -1))
       (n! k))))

#_
(defn nchoosek
  "Shermin idea for much cooler implementation, but at present about
   1/2 as fast...
   "
  [n k]
  (let [f (fn[n k]
            (reduce *' (take k (range n 1 -1))))]
    (/ (f n k) (f k k))))


(defn primes
  "RHickey paste.lisp.org with some fixes by jsa. Returns a list of
   all primes from 2 to n.  Wicked fast!"
  [n]
  (if (< n 2)
    ()
    (let [n (long n)]
      (let [root (long (Math/round (Math/floor (Math/sqrt n))))]
        (loop [i (long 3)
               a (int-array (inc n))
               result (list 2)]
          (if (> i n)
            (reverse result)
            (recur (+ i (long 2))
                   (if (<= i root)
                     (loop [arr a
                            inc (+ i i)
                            j (* i i)]
                       (if (> j n)
                         arr
                         (recur (do (aset arr j (long 1)) arr)
                                inc
                                (+ j inc))))
                     a)
                   (if (zero? (aget a i))
                     (conj result i)
                     result))))))))


(def +prime-set+ (atom (primes 1000)))

(defn prime-factors
  "Return the prime factorization of num as a seq of pairs [p n],
   where p is a prime and n is the number of times it is a factor.
   Ex: (prime-factors 510) => [[2 1] [3 1] [5 1] [17 1]]
  "
  [num]
  (if (< num 2)
    num
    (do
      (when (> num (last @+prime-set+))
        (swap! +prime-set+ (fn[_](primes (+ num (long (/ num 2)))))))
      (loop [ps (take-until #(> % num) @+prime-set+)
             factors []]
        (if (empty? ps)
          factors
          (let [pf (loop [n num
                          cnt 0]
                     (let [f (first ps)
                           [q r] (div n f)]
                       (if (not= 0 r)
                         (when (> cnt 0) [f cnt])
                     (recur q (inc cnt)))))]
            (recur (rest ps)
                   (if pf (conj factors pf) factors))))))))




;;; -----------------------------------------------------------------
;;; Simple vector stuff.  dot product, norm, distances, and such.  Is
;;; all this stuff in Incanter??


(defn dot [v1 v2]
  (double (reduce + 0.0 (map * v1 v2))))

(defn norm [v]
  (math/sqrt (dot v v)))

(defn vhat [v]
  (let [n (norm v)] (vec (map #(/ % n) v))))

(defn cos-vangle [v1 v2]
  (dot (vhat v1) (vhat v2)))

(defn vangle-dist [v1 v2]
  (math/abs (dec (cos-vangle v1 v2))))

(defn vecdist [v1 v2]
  (let [v (map - v1 v2)] (math/sqrt (dot v v))))

(defn vecmean
  ([v1 v2]
     (map #(/ (+ %1 %2) 2) v1 v2))
  ([vs]
     (let [d (count vs)]
       (map (fn[xs] (/ (apply + xs) d))
            (transpose vs)))))




;;; -----------------------------------------------------------------
;;; Various ancillary string stuff.  Again, eventually these pieces
;;; should be refactored into utils.* name spaces and files

;;; These are a bit duff.  Something like this should really be in
;;; clojure.string, clojure.contrib.string or clojure.str-utils
;;;
(defn string-less?
  "Case insensitve string comparison.  Usable as a sort comparator"
  [l r]
  (neg? (.compareToIgnoreCase l r)))

(defn string-greater?
  "Case insensitve string comparison.  Usable as a sort comparator"
  [l r]
  (pos? (.compareToIgnoreCase l r)))

(defn string-equal?
  "Case insensitve string comparison.  Usable as a sort comparator"
  [l r]
  (zero? (.compareToIgnoreCase l r)))


(defn intstg?
  "Test and return whether S is a string of only digits 0-9.  If so,
  return generalized boolean else return nil/false"
  [s]
  (let [hit (re-find #"[0-9]+" s)]
    (and hit (= hit s))))


(defn partition-stg
  "Returns a sequence of strings of n chars each, at offsets step
  apart. If step is not supplied, defaults to n, i.e. the partitions
  do not overlap. If a pad collection (a string or other collection of
  chars) is supplied, use its characters as necessary to complete last
  partition upto n characters. In case there are not enough padding
  chars, return a partition with less than n characters."
  ([n stg]
     (partition-stg n n stg))
  ([n step stg]
     (partition-stg n step "" stg))
  ([n step pad stg]
     (let [pad (if (string? pad) pad (str/join "" pad))
           stg (str stg pad)]
       (loop [s stg
              sv []]
         (if (= s "")
           sv
           (recur (str/drop step s)
                  (conj sv (str/take n s))))))))




;;; ----------------------------------------------------------------
;;; Slightly abstracted things from Java that are often used from/in
;;; many contexts...

(defn sys-property [prop-name]
  "Return the System property with name PROP-NAME (a string)"
  (System/getProperty prop-name))

(defn sys-properties []
  "Return the set of System properties as a map"
  (System/getProperties))

(defn classpath []
  (str/split #":" (sys-property "java.class.path")))

(defn getenv
  "Return the value of the environment variable EV (a string).
   Optionally, (no argument variant) return a map of all the
   environment variables.
  "
  ([ev] (System/getenv ev))
  ([] (System/getenv)))




;;; ----------------------------------------------------------------
;;; Some simple text file processing utils.


(defn reduce-file [fname func acc]
  "Reduce text file denoted by FNAME (filespec/File obj) using function
   FUNC per line and reduction seed accumulator ACC.  FUNC is a function
   of two arguments: first is the current line from file and second is
   the current value of ACC"
  (reduce func acc (io/read-lines (io/file-str fname))))


(defn process-line [acc line]
  (reduce #(assoc %1 %2 (inc (get %1 %2 0))) acc (str/split #"," line)))


(defmacro do-text-file [[in & options] & body]
  `(doseq [~'$line (io/read-lines (io/file-str ~in))]
     (do ~@body)))


(defmacro do-text-to-text [[in out] & body]
  `(io/with-out-writer (io/file-str ~out)
     (doseq [~'$line (io/read-lines (io/file-str ~in))]
       (let [result# (do ~@body)]
         (when result#
           (println result#))))))




;;; ----------------------------------------------------------------
;;; Definition macros and helpers for providing proper keyword
;;; arguments and maintaining doc string, user meta data, and special
;;; pre and post condition meta data.  Actually, this is just curried
;;; into the pre and post processing of the body.  This should really
;;; be resubmitted to clojure.contrib.def as the defnk and helper
;;; there do not account for pre and post conditions.


(defn name-with-attrs [name macro-args]
  "To be used in macro definitions.
   Handles optional docstrings and attribute maps for a name to be defined
   in a list of macro arguments. Also handles pre/post conditions. If the
   first macro argument is a string, it is added as a docstring to name and
   removed from the macro argument list. If afterwards the first macro argument
   is a map, its entries are added to the name's metadata map and the map is
   removed from the macro argument list. If the first form past the arg list
   is a map with :pre and/or :post, the map is removed and the pre and post
   vectors of forms are separated out.  The return value is a vector containing
   the name with its extended metadata map, the args form, followed by the pre
   and post forms (if any) and lastly, the body forms."
  [name macro-args]
  (let [[docstring macro-args] (if (string? (first macro-args))
                                 [(first macro-args) (next macro-args)]
                                 [nil macro-args])
        [attr macro-args]      (if (map? (first macro-args))
                                 [(first macro-args) (next macro-args)]
                                 [{} macro-args])
        attr                   (if docstring
                                 (assoc attr :doc docstring)
                                 attr)
        attr                   (if (meta name)
                                 (conj (meta name) attr)
                                 attr)
        [pre-post args body]   (if (and (map? (second macro-args))
                                        (or ((second macro-args) :pre)
                                            ((second macro-args) :post)))
                                 [(second macro-args)
                                  (first macro-args)
                                  (drop 2 macro-args)]
                                 [nil (first macro-args) (drop 1 macro-args)])
        [pre post]             (if pre-post
                                 [(pre-post :pre) (pre-post :post)]
                                 [nil nil])]
    [(with-meta name attr) args pre post body]))


(defmacro defnk [name & args-body]
  "Same as DEFN, but supports keyword style arguments.  Adapted and modified
   from Rich Hickey's original on google groups Clojure group, to support doc
   strings, meta data, and pre/post conditions."
  (let [[sym args pre post body] (name-with-attrs name args-body)
        pos-keys (split-with (complement keyword?) args)
        ps (pos-keys 0)
        ks (apply array-map (pos-keys 1))
        gkeys (gensym "gkeys__")
        letk (fn [ke]
               (let [k (key ke)
                     ;; The next oddity is due to some weird bug in
                     ;; Clojure 1.2 - for some reason in this context
                     ;; name returns nil despite k being a keyword!??!
                     kname (symbol (if (name k) (name k) (subs (str k) 1)))
                     v (val ke)]
                 `(~kname (if (contains? ~gkeys ~k) (~gkeys ~k) ~v))))]
    `(defn ~sym [~@ps & k#]
       (let [~gkeys (apply hash-map k#)
             ~@(apply concat (map letk ks))]
         ~@(if pre
             `((assert ~(conj (seq pre) 'and)))
             ())
         (let [res# (do ~@body)]
           ~@(if post
               `((assert ~(conj (map (fn [v] (replace `{% res#} v))
                                      (seq post))
                                'and)))
               ())
           res#)))))




;;; ----------------------------------------------------------------
;;; Some helpers for enhanced / simpler exception raising and
;;; handling.


;;; All of this is based on the support in slingshot lib and
;;; clojure.stacktrace.
;;;

(defmacro raise
  "Wraps throw+ for throwing a map object.  Ensures that (count args)
   is even, then places pairs into a map which is given to throw+
  "
  [& args]
  (let [m (into {} (map ensure-vec (partition 2 args)))]
    `(throw+ ~m)))

(defmacro handler-case
  "Hide some try+ details.  Form is the expression that is exception
   protected.  Cases are handler arms of the form [c e & body].  Where
   c is the predicate for the case arm (with out preceding 'catch'), e
   is the variable (symbol) to hold the exception object and body is
   the set of forms to execute.  Within body, the captured symbol
   contextMap is holds the exception context map.
  "
  [form & cases]
  `(try+
    ~form
    ~@(map (fn [[c# e# & body#]]
             `(catch ~c# ~e#
                (let [~'contextMap ~'&throw-context]
                  ~@body#)))
           cases)))


(defmacro with-handled
  "Wraps FORM in a slingshot try+ with catch arms for each condition
   in CONDITIONS. Each handle arm catches the condition C and prints a
   stack trace for it.  Hence, while this catches the conditions, it
   stops execution.  For catch and continue see CATCH-ALL
  "
  [form & conditions]
  `(try+
    ~form
    ~@(map (fn [c]
             `(catch [:type ~c] x#
                (print-stack-trace
                 (:throwable ~'&throw-context))))
           conditions)))


(defmacro with-ckd
  "FORMS wrapped by handlers and last ditch try/catch for standard exceptions.
   For conditions, the condition info meta data map is returned. For exceptions,
   the exception is converted to a string representation and that is returned."
  [& forms]
  `(try+
    (do ~@forms)
    (catch #(or (map? %) (set? %)) c#
      ~'&throw-context)
    (catch Exception e#
      (with-out-str
        (print e#)))))


(defmacro catch-all
  "FORMS wrapped by handlers and last ditch try/catch for standard exceptions.
   For conditions, the condition info meta data map is returned. For exceptions,
   the exception is converted to a string representation and that is returned."
  [& forms]
  `(with-ckd ~@forms))




;;; ----------------------------------------------------------------
;;; Some helpers for running external programs.  In particular,
;;; running them while ensuring they actually terminate normally.


(defn get-tool-path [toolset-type]
  (case toolset-type
        :ncbi
        (or (getenv "NCBI_BLAST")
            "/usr/local/ncbi/blast/bin/")
        :cmfinder
        (or (getenv "CMFINDER")
            "/usr/local/CMFinder/bin/")
        :infernal
        (or (getenv "INFERNAL")
            "/usr/local/Infernal/bin/")

        :cdhit
        (or (getenv "CDHIT")
            "/usr/local/cd-hit/")

        :bioperl "/usr/local/bin/"
        :mysql   "/usr/local/mysql/bin/"))


(defn assert-tools-exist [tool-vec]
  (doseq [pgm (vec tool-vec)]
    (when (not (fs/executable? pgm))
      (let [[_ path p] (re-matches #"(^/.*/)(.*)" pgm)]
        (raise :type :missing-program :path path :pgm p)))))


(defn runx
  "RUNs an eXternal program PROGRAM (a string naming the program executable),
   passing it the (typically, _command line_) arguments ARGS.  ARGS is either
   a single vector or list of the arguments, or a sequential list of the
   arguments given to the function as further arguments after the program."
  [program & args]
  (let [the-args (vec (if (and (= 1 (count args))
                               (sequential? (first args)))
                        (first args)
                        args))

        i (first (seq/positions  #(= :> %1) the-args))

        [the-args stdout-file] (if (not i)
                                 [the-args nil]
                                 [(vec (concat (take i the-args)
                                               (subvec the-args (+ 2 i))))
                                  (the-args (+ 1 i))])

        i (first (seq/positions #(= :?> %1) the-args))
        [the-args stderr-file] (if (not i)
                                 [the-args nil]
                                 [(concat (take i the-args)
                                          (subvec the-args (+ 2 i)))
                                  (the-args (+ 1 i))])

        write-out (fn[stdout]
                    (let [rdr (io/reader stdout)]
                      (io/with-out-writer (fs/fullpath stdout-file)
                        (doseq [l (line-seq rdr)] (println l)))))
        write-err (fn[stderr]
                    (let [rdr (io/reader stderr)]
                      (io/with-out-writer (fs/fullpath stderr-file)
                        (doseq [l (line-seq rdr)] (println l)))))

        the-args (if stdout-file (concat the-args [:out-fn write-out]) the-args)
        the-args (if stderr-file (concat the-args [:err-fn write-err]) the-args)

        result (apply sh program the-args)
        result (if stderr-file (assoc result :err stderr-file) result)
        result (if stdout-file (assoc result :out stdout-file) result)]

    (when (not= 0 (result :exit))
      (if stderr-file
        :err-in-result
        (raise :type :program-failed
               :exit (result :exit)
               :pgm program :err (result :err)
               :args the-args)))
    (if stdout-file
      result
      (result :out))))
