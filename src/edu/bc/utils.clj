(ns edu.bc.utils

  "General utility functions and macros.  Basically these resources are
   fairly general and intended to be usable on most any project"

  (:require [clojure.contrib.math :as math]
            [clojure.contrib.combinatorics :as comb]
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.contrib.properties :as prop]
            [clojure.xml :as xml]
            [edu.bc.fs :as fs])

  (:use [clojure.contrib.condition
         :only [raise handler-case *condition*
                print-stack-trace stack-trace-info]]
        [clj-shell.shell
         :only (sh)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)])

  (:import (java.util Date Calendar Locale)
           java.lang.Thread
           (java.text SimpleDateFormat)))



;;; (compile 'edu.bc.utils)


;;; -----------------------------------------------------------------
;;; Various ancillary math/arithmetic stuff.

(defn div
  "Integer division.  Return [q r], such that floor(n / d) * q + r = q"
  [n d]
  (let [q (math/floor (/ n d))]
    [q (rem n d)]))


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
  (if (or (= n 0) (= n 1))
    1
    (reduce * (range 2 (inc n)))))

(defn n-k!
  "For positive integers N and K, compute (n-k)!"
  [n k]
  {:pre [(integer? n) (integer? k) (> n -1) (<= 0 k n)]}
  (if (= n k)
    1
    (n! (- n k))))

(defn nCk
  "For positive integers N and K, compute N choose K (binomial
   coefficient): n!/((n-k)!k!)"
  [n k]
  {:pre [(integer? n) (integer? k) (> n -1) (<= 0 k n)]}
  (/ (reduce * (range n (- n k) -1))
     (n! k)))


(defn shannon-entropy
  "Returns the Shannon entropy of a sequence: -sum(* pi (log pi)),
   where i ranges over the unique elements of S and pi is the
   probability of i in S: (freq i s)/(count s)"
  [s]
  (let [fs (frequencies s)
        cnt (float (reduce (fn[x [k v]] (+ x v)) 0 fs))
        H (reduce
           (fn[H [c v]]
             (let [p (float (/ v cnt))]
               (float (+ H (float (* p (float (log2 p))))))))
           (float 0.0) fs)]
    (float (- H))))




(defn sum
  "Return the sum of the numbers in COLL.  If COLL is a map, return
  the sum of the (presumed all numbers) in (vals coll)"
  [coll]
  (let [v (if (map? coll) (vals coll) coll)]
    (apply + v)))


(defn primes
  "RHickey paste.lisp.org with some fixes by jsa. Returns a list of
   all primes from 2 to n.  Wicked fast!" [n]
  (if (< n 2)
    ()
    (let [n (int n)]
      (let [root (int (Math/round (Math/floor (Math/sqrt n))))]
        (loop [i (int 3)
               a (int-array (inc n))
               result (list 2)]
          (if (> i n)
            (reverse result)
            (recur (+ i (int 2))
                   (if (<= i root)
                     (loop [arr a
                            inc (+ i i)
                            j (* i i)]
                       (if (> j n)
                         arr
                         (recur (do (aset arr j (int 1)) arr)
                                inc
                                (+ j inc))))
                     a)
                   (if (zero? (aget a i))
                     (conj result i)
                     result))))))))


(defn prime-factors [num]
  (loop [ps (primes (int (math/sqrt num)))
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
               (if pf (conj factors pf) factors))))))




;;; -----------------------------------------------------------------
;;; Miscellaneous changes/additions/"fixes" to various operations that
;;; are either already in Clojure or in bits of contribs or probably
;;; _will_ be in future (at which point these can be retired...)


(defn drop-until [pred coll] (drop-while (complement pred) coll))
(defn take-until [pred coll] (take-while (complement pred) coll))


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


(defn ensure-vec [x]
  (cond
   (vector? x) x
   (seq? x) (vec x)
   (map? x) (vec x)
   true [x]))


(defn in [e coll]
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


(defn random-subset [s cnt]
  (let [s (vec (set s))]
    (cond
     (<= cnt 0) #{}
     (<= (count s) cnt) (set s)
     :else
     (loop [ss #{(rand-nth s)}]
       (if (= (count ss) cnt)
         ss
         (recur (conj ss (rand-nth s))))))))

(defn combins [k coll]
  (lazy-seq (map vec (comb/combinations coll k))))

(defn choose-k [k coll]
  (combins k coll))





(def ^{:private true} *uid* (atom (.getTime (Date.))))

(defn gen-uid []
  (swap! *uid* inc))

(defn gen-kwuid []
  (keyword (str (gen-uid))))



(defn sleep [msecs]
  (Thread/sleep msecs))

(defn str-date
  ([] (str-date (Date.) "yyyy-MM-dd HH:mm:ss"))
  ([fm] (str-date (Date.) fm))
  ([d fm] (.format (SimpleDateFormat. fm) d)))




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

(defn levenshtein [s t]
  (cond
   (= s t "") 0
   (= 0 (count s)) (count t)
   (= 0 (count t)) (count s)
   :else
   (peek (reduce
          (fn [prev-row s-elem] (new-row prev-row s-elem t))
          (range (inc (count t)))
          s))))




;;; Various frequency counts, probabilities, similarity coefficients,
;;; corresponding difference fns and ngram operations

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

(defn freqs-probs [n coll]
  (let [freqs (freqn n coll)
        sz (sum (vals freqs))
        probs (reduce (fn[m [k v]]
                        (assoc m k (float (/ v sz))))
                      {} freqs)]
    [freqs probs sz]))


(defn cc-freqn [n colls]
  (reduce
   (fn[m coll]
     (merge-with #(+ %1 %2) m (freqn n coll)))
   {} colls))

(defn cc-freqs-probs [n colls]
  (let [freqs (cc-freqn n colls)
        sz (sum (vals freqs))
        probs (reduce (fn[m [k v]]
                        (assoc m k (float (/ v sz))))
                      {} freqs)]
    [freqs probs sz]))


(defn alphabet [coll & {n :n :or {n 1}}]
  (cond
   (map? coll) (keys coll)
   (set? coll) coll
   :else
   (keys (freqn n coll))))

(defn combins-freqn [n coll]
  (reduce (fn[m x] (assoc m x (inc (get m x 0))))
          {} (combins n coll)))

(defn choose-k-freqn [n coll] )




(defn probs [n coll]
  (second (freqs-probs n coll)))


(defn log-odds [frq1 frq2]
  (ln (/ frq1 frq2)))

(defn lod-score [qij pi pj]
  (ln (/ qij (* pi pj))))

(defn raw-lod-score [qij pi pj & {scaling :scaling :or {scaling 1.0}}]
  (if (= scaling 1.0)
    (lod-score qij pi pj)
    (int (/ (lod-score qij pi pj) scaling))))

(defn bg-freq [n seqs] )

(defn bg-freqs-probs [] )








(defn diff-fn
  "Return the function that is 1 - F applied to its args: (1 - (apply f args)).

   Ex: (let [dice-diff (diff-fn dice-coeff)
             ...]
         (dice-diff some-set1 some-set2))"
  [f] (fn [& args] (- 1 (apply f args))))


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
        c1 (freq-sum (set/intersection freq-s1 freq-s2))
        c2 (freq-sum (set/union freq-s1 freq-s2))]
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
  [s1 s2 & {uc? :uc? n :n sfn :sfn ngfn :ngfn
            :or {n 2 uc? false sfn dice-coeff ngfn word-letter-pairs}}]
  (let [s1 (ngfn n (if uc? (str/upper-case s1) s1))
        s2 (ngfn n (if uc? (str/upper-case s2) s2))]
    (sfn s1 s2)))

;;;(float (ngram-compare "FRANCE" "french"))
;;;(float (ngram-compare "FRANCE" "REPUBLIC OF FRANCE"))
;;;(float (ngram-compare "FRENCH REPUBLIC" "republic of france"))
;;;(float (ngram-compare
;;;        "TATATTTGGAGTTATACTATGTCTCTAAGCACTGAAGCAAA"
;;;        "TATATATTTTGGAGATGCACAT"))




;;; -----------------------------------------------------------------
;;; Simple vector stuff.  dot product, norm, distances, and such.


(defn dot [v1 v2]
  (reduce #(+ %1 %2) 0 (map #(* %1 %2) v1 v2)))

(defn norm [v]
  (math/sqrt (dot v v)))

(defn vhat [v]
  (let [n (norm v)] (vec (map #(/ % n) v))))

(defn cos-vangle [v1 v2]
  (dot (vhat v1) (vhat v2)))

(defn vangle-dist [v1 v2]
  (math/abs (dec (cos-vangle v1 v2))))

(defn vecdist [v1 v2]
  (let [v (vec (map #(- %1 %2) v1 v2))] (dot v v)))


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

(defn getenv [ev]
  "Return the value of the environment variable EV (a string)"
  (System/getenv ev))




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
;;; Some helpers for running external programs.  In particular,
;;; running them while ensuring they actually terminate normally.


;;; The "odd" explicit use of list is the result of a toxic
;;; interaction between syntax-quote not generating lists (but cons
;;; cells) and handler-case not taking this into account.  The bug is
;;; clearly in handler-case as it should account for this
;;;
(defmacro with-handled
  "Wraps FORM in a handler-case with handle case arms for each condition in
   CONDITIONS. Each handle arm catches the condition C and prints a stack trace
   for it.  Hence, while this catches the conditions, it stops execution.
   For catch and continue see CATCH-ALL"
  [form & conditions]
  `(handler-case :type
     ~form
     ~@(map (fn [c] (list 'handle c
                           `(stack-trace-info *condition*)))
            conditions)))


(defmacro with-ckd
  "FORMS wrapped by handlers and last ditch try/catch for standard exceptions.
   For conditions, the condition info meta data map is returned. For exceptions,
   the exception is converted to a string representation and that is returned."
  [& forms]
  `(try
     (do ~@forms)
     (catch clojure.contrib.condition.Condition c#
       (meta c#))
     (catch Exception e#
       (with-out-str
         (print e#)))))


(defmacro catch-all
  "FORMS wrapped by handlers and last ditch try/catch for standard exceptions.
   For conditions, the condition info meta data map is returned. For exceptions,
   the exception is converted to a string representation and that is returned."
  [& forms]
  `(with-ckd ~@forms))



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
