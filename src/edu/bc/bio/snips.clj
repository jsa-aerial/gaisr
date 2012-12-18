
(in-ns 'edu.bc.bio.gaisr.post-db-csv)


(defn pval2
  "Reads in a sto file and n = an integer. The sto file is used to
   generate a Clustal W type file with the same file name with an .aln
   extension instead of .sto.  Produces .aln and generates n random
   shuffled alignments and then calcuates p value by taking the
   fraction sampled alignment MI > true MI. returns a position-pair
   and the pvalue"

  [stoin n]
  (let [aln (sto->aln stoin (str (subs stoin 0 (- (count stoin) 3)) "aln"))
        p (profile (read-sto stoin))
        mi (mutual_info p)
        rand_mi (apply concat (pmap (fn [randa]
                                      (map #(mutual_info (profile %)) randa))
                                    (partition-all (/ n 4) (rand_aln aln n))))]
    (for [k (keys mi)]
      [k (double (/ (count
                     (filter (fn [x]
                               (let [sample-mi (get x k)]
                                 (> sample-mi (get mi k))))
                             rand_mi))
                    n))])))





;;; Sto to ALN stuff
;;; NOTE: This is in edu.bc.bio.gaisr.post-db-csv NAMESPACE
;;;
;;; Create them
(let [base "/data2/Bio/Training/MoStos2"
      dirs (fs/directory-files base "")
      pos-stos (map #(fs/directory-files % "pos.sto") dirs)
      neg-stos (map #(fs/directory-files % "neg.sto") dirs)
      dir-stos (map #(concat %1 %2) pos-stos neg-stos)]
  (doseq [dstos dir-stos]
    (doseq [ds dstos]
      (when (not (fs/empty? (fs/replace-type ds ".fna")))
        (sto->aln ds (fs/replace-type ds ".aln"))))))

;;; Count them all up
(let [base "/data2/Bio/Training/MoStos"
      dirs (fs/directory-files base "")
      pos-alns (flatten (map #(fs/directory-files % "pos.aln") dirs))
      neg-alns (flatten (map #(fs/directory-files % "neg.aln") dirs))]
  [(count (concat pos-alns neg-alns)) (count pos-alns) (count neg-alns)])
==> [145 61 84]

;;; Grab a couple +'s and -'s to look at
(let [base "/data2/Bio/Training/MoStos"
      dirs (fs/directory-files base "")
      pos-alns (flatten (map #(fs/directory-files % "pos.aln") dirs))
      neg-alns (flatten (map #(fs/directory-files % "neg.aln") dirs))]
  [(take 2 pos-alns) (take 2 neg-alns)])

;;; Split up for xvalidation (counting the pos and negs in partitions)
(let [folds 5
      base "/data2/Bio/Training/MoStos"
      dirs (fs/directory-files base "")
      pos-alns (flatten (map #(fs/directory-files % "pos.aln") dirs))
      neg-alns (flatten (map #(fs/directory-files % "neg.aln") dirs))
      both (concat pos-alns neg-alns)
      size (int (/ (count both) folds))
      sets (partition size size [] (shuffle both))
      sets (map (fn[x]
                  (let [[p n] (seq/separate #(re-find #"pos" %) x)]
                    [(count p) (count n)])) sets)]
  (prn :*** size sets)
  (doseq [r (seq/rotations sets)]
    (prn (first r) (rest r))))

(def folds
     (let [folds 5
           base "/data2/Bio/Training/MoStos"
           dirs (fs/directory-files base "")
           pos-alns (flatten (map #(fs/directory-files % "pos.aln") dirs))
           neg-alns (flatten (map #(fs/directory-files % "neg.aln") dirs))
           both (concat pos-alns neg-alns)
           size (int (/ (count both) folds))
           sets (partition size size [] (shuffle both))
           sets2 (map (fn[x]
                        (let [[p n] (seq/separate #(re-find #"pos" %) x)]
                          [(count p) (count n)])) sets)]
       (prn :*** size)
       (map #(do [(first %) (rest %)])
            (seq/rotations sets))))


(def folds
     (let [folds 10
           base "/data2/Bio/Training/MoStos2"
           dirs (fs/directory-files base "")
           pos-alns (flatten (map #(fs/directory-files % "pos.aln") dirs))
           neg-alns (flatten (map #(fs/directory-files % "neg.aln") dirs))
           both (concat pos-alns neg-alns)
           [both-firms both-proteo] (seq/separate #(re-find #"firm" %) both)
           testset both
           size (int (/ (count testset) folds))
           sets (partition size size [] (shuffle testset))
           sets2 (map (fn[x]
                        (let [[p n] (seq/separate #(re-find #"pos" %) x)]
                          [(count p) (count n)])) sets)]
       (prn :*** size)
       (map #(do [(first %) (rest %)])
            (seq/rotations sets))))


(time
 (doseq [fold folds]
   (prn (catch-all (let [resp (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel"
                                    (concat ["-n" "/home/jsa/TMP/km.dat"]
                                            (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                                          (seq/separate #(re-find #"pos" %) (flatten (second fold)))))))
                         tra (runx "/usr/local/libsvm/svm-train" "-t" "4" "-b" "1" "/home/jsa/TMP/km.dat" "/home/jsa/TMP/km.model")
                         predx (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel"
                                     (concat ["-n" "/home/jsa/TMP/x.dat"]
                                             (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                                           (seq/separate #(re-find #"pos" %) (flatten (second fold)))))
                                             ["--test"]
                                             (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                                           (seq/separate #(re-find #"pos" %) (map #(fs/replace-type % ".fna") (first fold)))))))
                         preds (runx "/usr/local/libsvm/svm-predict" "-b" "1" "/home/jsa/TMP/x.dat" "/home/jsa/TMP/km.model" "/home/jsa/TMP/x.output")]
                     preds)))))

(def tst-90 (ffirst folds))
(def tr-84 (second (nth folds 9)))




;;; String BPLA kernel training and SVM prediction runs
;;;
(catch-all
 (let [resp (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel"
                  (concat ["-n" "/home/jsa/TMP/km.dat"]
                          (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                        (seq/separate #(re-find #"pos" %) (flatten (second (first folds))))))))
       tra (runx "/usr/local/libsvm/svm-train" "-t" "4" "-b" "1" "/home/jsa/TMP/km.dat" "/home/jsa/TMP/km.model")
       predx (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel"
                   (concat ["-n" "/home/jsa/TMP/x.dat"]
                           (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                         (seq/separate #(re-find #"pos" %) (flatten (second (first folds))))))
                           ["--test"]
                           (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                         (seq/separate #(re-find #"pos" %) (map #(fs/replace-type % ".fna") (ffirst folds)))))))
       preds (runx "/usr/local/libsvm/svm-predict" "-b" "1" "/home/jsa/TMP/x.dat" "/home/jsa/TMP/km.model" "/home/jsa/TMP/x.output")]
   preds))


(catch-all
 (let [resp (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel" (concat ["-n" "/home/jsa/TMP/km.dat"] (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s)) (seq/separate #(re-find #"pos" %) (ffirst folds))))))
       tra (runx "/usr/local/libsvm/svm-train" "-t" "4" "-b" "1" "/home/jsa/TMP/km.dat" "/home/jsa/TMP/km.model")
       predx (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel" (concat ["-n" "/home/jsa/TMP/x.dat"] (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s)) (seq/separate #(re-find #"pos" %) (ffirst folds)))) ["--test" "+1" "/data2/Bio/Training/MoStos/ECS4/ECS4_071211.sto.proteobacteria-genus-pos.fna" "-1" "/data2/Bio/Training/MoStos/ECS4/ECS4_071211.sto.firmicutes-genus-neg.fna"]))
       preds (runx "/usr/local/libsvm/svm-predict" "-b" "1" "/home/jsa/TMP/x.dat" "/home/jsa/TMP/km.model" "/home/jsa/TMP/x.output")]
   preds))


;;; Full 5 fold validation run
(time
 (doseq [fold folds]
   (prn (catch-all (let [resp (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel"
                                    (concat ["-n" "/home/jsa/TMP/km.dat"]
                                            (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                                          (seq/separate #(re-find #"pos" %) (flatten (second fold)))))))
                         tra (runx "/usr/local/libsvm/svm-train" "-t" "4" "-b" "1" "/home/jsa/TMP/km.dat" "/home/jsa/TMP/km.model")
                         predx (runx "/home/jsa/Clojure/mlab/AuxSoftWare/bpla_kernel-1.0/bpla_kernel/bpla_kernel"
                                     (concat ["-n" "/home/jsa/TMP/x.dat"]
                                             (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                                           (seq/separate #(re-find #"pos" %) (flatten (second fold)))))
                                             ["--test"]
                                             (flatten (map (fn[s] (map #(if (re-find #"pos" %) ["+1" %] ["-1" %]) s))
                                                           (seq/separate #(re-find #"pos" %) (map #(fs/replace-type % ".fna") (first fold)))))))
                         preds (runx "/usr/local/libsvm/svm-predict" "-b" "1" "/home/jsa/TMP/x.dat" "/home/jsa/TMP/km.model" "/home/jsa/TMP/x.output")]
                     preds)))))





;;; New metag1 cmsearches...
(def metag1-fasta-files
     (sort (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")))

(def new-metag1-cms
     (map #(first (fs/directory-files
                   (fs/join "/data2/Bio/MetaG1/MoStos" %) ".cm"))
          (str/split #" " "RF00050 RF00059 RF00168 RF00504")))

(cms&hitfnas->cmsearch-out
 new-metag1-cms metag1-fasta-files
 :cmpar 4 :eval 1.0)



;;; insert an origin line for stos without them...
;;;
(map #(join-sto-fasta-file
       % % :origin "#=GF AU Infernal 1.0.2")
     (sort (fs/directory-files "/data2/Bio/RFAM" ".sto")))


(let [base "/data2/Bio/Training/MoStos"
      RFAM  "/data2/Bio/RFAM"
      dir-f-pairs
      (keep (fn[f]
              (let [mosto-dir (fs/join base (first (str/split #"-" f)))]
                (when (not (fs/directory? mosto-dir)) [mosto-dir f])))
            (map fs/basename (sort (fs/directory-files RFAM ".sto"))))]
  (doseq [[d f] dir-f-pairs]
    (fs/mkdir d)
    (fs/copy (fs/join RFAM f) (fs/join d f))))

(def *cur-search*
     (future (catch-all
              (doall (cms&hitfnas->cmsearch-out
                      "/data2/Bio/Training/MoStos/RF01054/RF01054.sto.cm"
                      (directory-files "/data2/Bio/Training/FastaFiles" ".fna")
                      :eval 1.0)))))

(def *cur-search*
     (future (doall (cms&hitfnas->cmsearch-out
                     @*cur-run*
                     (sort (fs/directory-files
                            "/data2/Bio/Training/FastaFiles" ".fna"))
                     :eval 1.0))))


(pds/gen-aligned-training-sets
 "/data2/Bio/Training/MoStos/ECS4" :ev-value 0.00001)

(def *dirs-needing-training-sets*
     (let [base "/data2/Bio/Training/MoStos2"
           dirs (sort (filter fs/directory? (fs/directory-files base "")))
           names (map fs/basename dirs)]
       (filter (fn [d]
                 (and (not-empty (fs/directory-files d ".cmsearch.out"))
                      (reduce #(and %1 (not (fs/empty? %2)))
                              true (fs/directory-files d ".cmsearch.out"))
                      (empty? (fs/directory-files d ".cmsearch.csv"))))
               dirs)))


(def *dirs-needing-cms*
     (let [base "/data2/Bio/Training/MoStos"
           dirs (sort (filter fs/directory? (fs/directory-files base "")))
           names (map fs/basename dirs)]
       (filter (fn [d] (and (not-empty (fs/directory-files d ".sto"))
                            (empty? (fs/directory-files d ".cm"))))
               dirs)))

(def *cur-run*
     (future
      (mostos->calibrated-cms
       (map #(first (fs/directory-files % ".sto"))
            *dirs-needing-cms*))))




(def jdkbin (fs/listdir "/opt/Java/java/bin/"))
(def altbin (fs/listdir "/etc/alternatives/"))

(io/with-out-writer (fs/fullpath "~/fix-java" )
  (println "#!/bin/tcsh")
  (println)
  (doseq [x (map #(str "rm -f " %)
                 (sort (set/intersection (set jdkbin) (set altbin))))]
    (println x))
  (doseq [x (let [base "/opt/Java/java/bin/"]
              (map #(str "ln -s " base % " " %)
                   (sort (set/intersection (set jdkbin) (set altbin)))))]
    (println x)))



(def *jsa* nil)

(def *jsa*
  (let [stmt "select name,description from bioentry"]
    (sql/with-connection mysql-ds
      (sql/with-query-results qresults [stmt]
        (doall qresults)))))




(doseq [i (range 1 133)]
  (let [stmt (str "select t.name,sfqv.* from term as t, seqfeature_qualifier_value as sfqv where t.term_id=sfqv.term_id and sfqv.term_id=" i " limit 10;")]
    (sql/with-connection mysql-ds
      (sql/with-query-results qresults [stmt]
        (println "\n" stmt ":")
        (doseq [m (doall qresults)]
          (println i (m :name) "==>" (m :seqfeature_id) "," (m :value)))))))



(first *jsa*)
(count *jsa*)

(doseq [x (take 2 *jsa*)]
  (println (:name x) "\t" (:description x)))

(io/with-out-writer "/data2/Bio/ECRibLeaders/nc-names-only.txt"
  (doseq [x (sort-by :name *jsa*)]
    (println (:name x))))

(io/with-out-writer "/data2/Bio/ECRibLeaders/nc-names-descs.txt"
  (doseq [x (sort-by :name *jsa*)]
    (println (:name x) "\t" (:description x))))


(io/with-out-writer "/data2/Bio/ECRibLeaders/all-names-only.txt"
  (doseq [x (sort-by :name *jsa*)]
    (println (:name x))))

(io/with-out-writer "/data2/Bio/ECRibLeaders/all-names-descs.txt"
  (doseq [x (sort-by :name *jsa*)]
    (println (:name x) "\t" (:description x))))






(defn make-map [inseq]
  (reduce (fn[m [k v]]
            (assoc m k (conj (get m k []) v)))
          {} inseq))


(make-map [[:a 1] [:b 2] [:c 3] [:d 4] [:d 5] [:d 7]])


(defn process-filex [file]
  (make-map
   (map #(str/split #" " %)
        (io/read-lines file))))

(map #(str/split #" " %)
     (io/read-lines "/home/jsa/Bio/Blasting/CMFtest/Clostridium-203/jsa.txt"))

(process-filex "/home/jsa/Bio/Blasting/CMFtest/Clostridium-203/jsa.txt")


(defn operon-filter [entries & {par :par :or {par 10}}]
  (let [q (math/floor (/ (count entries) par))
        entsets (partition q q [] entries)
        leftover (last entsets)
        entsets (conj (drop 1 (butlast entsets))
                      (set/union (first entsets) leftover))]
    (reduce
     (fn[v subv]
       (set/union v subv))
     []
     (pmap (fn[entset]
             (doall
              (map (fn[x]
                     (let [[name loc] (str/split #" " x)
                           [s e] (map #(Integer. %) (str/split #"-" loc))
                           strand (if (< e s) -1 1)
                           [ns ne] (get-region name strand s)
                           diff (math/abs (- ns ne))]
                       [name strand [ns ne] diff]))
                   entset)))
           entsets))))


(def *test-set*
     (-> "/home/jsa/Bio/Blasting/b-subtilus-Ls.fna.blast"
         slurp
         ((fn[x](str/split #"\n" x)))
         (nlsq-tuples-from-hits "/home/jsa/Bio/Blasting/b-subtilus-Ls.hitfna")))

(def *test-map*
     (reduce (fn[m [n r sq]] (assoc m n (conj (get m n []) [n r sq])))
             {} *test-set*))


(def phylo-clusters (group-by-taxon *test-map*))

(def testsqs (let [[tx cnt tuples] (first phylo-clusters)] tuples))


(cd-hit-est testsqs "/home/jsa/Bio/Blasting/clj-cdhit.fna")
(def cdhitsqs
     (map (fn[[n sq]] [(subs n 1) sq])
          (partition
           2 (str/split #"\n" (slurp "/home/jsa/Bio/Blasting/clj-cdhit.fna")))))

(time (get-keeper-set testsqs :cmpfn ngram-closeness))
(map (fn[[n c sq]] [n sq]) *1)
(set/intersection (set *1) (set cdhitsqs))
(set/difference (set *1) (set cdhitsqs))



(defn get-keeper-set-orig [nlsq-tuples compare-set keepers]
  (if (empty? nlsq-tuples)
    keepers
    (let [q (first nlsq-tuples)
          test-set (remove #(= q %) compare-set)
          same? (seq-same? q test-set :cmpfn *cmpfn*)]
      (recur
       (drop 1 nlsq-tuples)
       (if same? (remove #(= q %) compare-set) compare-set)
       (if same? keepers (conj keepers q))))))




(def x (group-by second
                 (map (fn[m]
                        [(m :name)
                         (nth (reverse (str/split #", " (m :ancestors))) 3)])
                      (mlab/names->tax *test-map*))))

(reduce (fn [m [n tx]]
          (assoc m tx (set/union (get m tx #{}) (set (*test-map* n)))))
        {} (val (nth (seq x) 5)))


(map (fn[tuples]
       [(count tuples)
        (count (get-keeper-set tuples (set tuples ) #{}))])
     (list (vec (*1 "Micromonosporaceae"))))



(map (fn[tuples]
       [(count tuples)
        (count (get-keeper-set tuples (set tuples ) #{}))])
     (list (vec (*1 "Polynucleobacter necessarius"))))



(defn array? [x] (-> x class .isArray))
(defn see [x] (if (array? x) (map see x) x))

(see (make-array Float/TYPE 100))
(see (into-array Float/TYPE [1 2 3 4 5 6 7 8 9 10]))

(def sfa (make-array Float/TYPE 1000000))


(defmacro %aget
  ([hint array idx]
    `(aget ~(vary-meta array assoc :tag hint) (int ~idx)))
  ([hint array idx & idxs]
    `(let [a# (aget ~(vary-meta array assoc :tag 'objects) (int ~idx))]
       (%aget ~hint a# ~@idxs))))

(defmacro %aset [hint array & idxsv]
  (let [hints '{floats float doubles double ints int}
        [v idx & sxdi] (reverse idxsv)
        idxs (reverse sxdi)
        v (if-let [h (hints hint)] (list h v) v)
        nested-array (if (seq idxs)
                       `(%aget ~'objects ~array ~@idxs)
                        array)
        a-sym (with-meta (gensym "a") {:tag hint})]
      `(let [~a-sym ~nested-array]
         (aset ~a-sym (int ~idx) ~v))))


(%aget floats sfa 9)
(%aset floats sfa 9 Math/E)
(%aset floats sfa 9 Math/PI)
(%aget floats sfa (int 9))

(defn #^Float xx [#^Float x #^Float y]
  (+ (float x) (float y)))

(defn yy [#^Float x #^Float y]
  (+ (float x) (float y)))

(defn foo [x y]
  (+ (float 4.5) (yy x y)))



(defn #^Float dopoly [#^Float x]
  (let [x (float x)
        cnt (int 100)
        pol (make-array Float/TYPE 100)]
    (loop [j (int 0)
           mu (float 10.0)]
        (when (< j cnt)
          (%aset floats pol j (/ (+ mu (float 2.0)) (float 2.0)))
          (recur (unchecked-inc j)
                 (%aget floats pol j))))
    (loop [j (int 0)
           su (float 0.0)]
        (if (< j cnt)
          (recur (unchecked-inc j)
                 (+ (%aget floats pol j) (* su x)))
          su))))


(defn eval-pol [#^Integer n #^Float x]
  (let [n (int n)
        x (float x)]
    (loop [i (int 0)
           pu (float 0.0)]
      (if (< i n)
        (recur (unchecked-inc i)
               (+ pu (float (dopoly x))))
        pu))))


(defn eval-pol-par [#^Integer n #^Float x & {par :par :or {par 10}}]
  (let [n (int n)
        q (math/floor (/ n par))
        r (rem n 10)
        chunks (conj (repeat (dec par) q) (+ q r))]
    (reduce #(+ (float %1) (float %2))
            (float 0.0)
            (pmap #(eval-pol % x) chunks))))



(map #(do [% (+ % 1000000 -1)])
     (take 10 (iterate #(+ % 1000000) 1)))



(defn div-3-5 [n] (or (= 0 (rem n 3)) (= 0 (rem n 5))))

(defn check-chunk [s e]
  (count (filter div-3-5 (range s (inc e)))))

(defn par-count []
  (reduce
   + 0 (pmap (fn[[x y]] (check-chunk x y))
             (map #(do [% (+ % 1000000 -1)])
                  (take 10 (iterate #(+ % 1000000) 1))))))




;;;-----------------------------------------------------------------------

(io/with-out-writer "/home/jsa/Bio/bigrun-all.csv"
  (println edu.bc.bio.gaisr.post-db-csv/+cmsearch-csv-header+)
  (doseq [l (->> (fs/directory-files
		  "/home/kaila/Bio/STOfiles/RNA-082712/BigRun/CSV"
		  ".cmsearch.csv")
		 (filter #(> (count (io/read-lines %)) 1))
		 (reduce (fn[v f] (->> f io/read-lines (drop 1) (concat v)))
			 []))] 
    (println l)))


;;; ----------------------------------------------------------------------

(def *bg-freqs-probs-all-rfam-alns-mon*
     (bg-freqs-probs 1 "/data2/Bio/RFAM" :ftypes [".sto"] :norm true :par 2))

(def *bg-freqs-probs-all-rfam-alns-di*
     (bg-freqs-probs 2 "/data2/Bio/RFAM" :ftypes [".sto"] :norm true :par 2))

(def *bg-freqs-probs-all-rfam-alns-col*
     (bg-freqs-probs 2 "/data2/Bio/RFAM"
                     :fsps-fn cc-combins-freqs-probs
                     :cols true :ftypes [".sto"] :norm true :par 2))

;;; Assert symmetric
(def *bg-freqs-probs-all-rfam-alns-col*
     [(coalesce-xy-yx
       (sort-by second > (first *bg-freqs-probs-all-rfam-alns-col*))
       (fn [x v]
         (if (not v) 0 (+ (val x) v))))
      (coalesce-xy-yx
       (sort-by second > (second *bg-freqs-probs-all-rfam-alns-col*))
       (fn [x v]
         (if (not v) 0 (+ (val x) v))))])

(def *bg-freqs-probs-all-rfam-alns-col-nonsym*
     (bg-freqs-probs 2 "/data2/Bio/RFAM"
                     :fsps-fn cc-combins-freqs-probs
                     :cols true :ftypes [".sto"] :norm true :par 2))


(sort-by second > (first *bg-freqs-probs-all-rfam-alns-mon*))
(sort-by second > (second *bg-freqs-probs-all-rfam-alns-mon*))

(sort-by second > (first *bg-freqs-probs-all-rfam-alns-di*))
(sort-by second > (second *bg-freqs-probs-all-rfam-alns-di*))

(sort-by second > (first *bg-freqs-probs-all-rfam-alns-col*))
(sort-by second > (second *bg-freqs-probs-all-rfam-alns-col*))

(defn print-dist [dist]
  (doseq [[k v] (sort-by second > dist)]
    (print " ") (pr k) (print " ") (pr v) (println \,)))








Corresponding probabilities

{"GC" 0.08437428795109882,
 "CG" 0.08071853660437028,
 "TT" 0.0696883223046914,
 "AA" 0.06967938106879766,
 "CC" 0.06552183985818097,
 "GG" 0.06547925910697425,
 "AT" 0.06386396929374653,
 "CA" 0.06223383629426149,
 "TG" 0.06220668668188498,
 "TC" 0.06140405325405737,
 "GA" 0.06139930008570952,
 "AG" 0.05361838101285567,
 "CT" 0.05361308751758834,
 "AC" 0.05078720100758968,
 "GT" 0.05077003569698625,
 "TA" 0.04463641920557993,
 "NN" 2.530566667579171E-6,
 "YT" 1.480025909893127E-7,
 "NA" 1.46934751227485E-7,
 "AN" 1.443719357990986E-7,
 "TN" 1.418091203707123E-7,
 "NT" 1.336935381808221E-7,
 "NG" 1.025126171354547E-7,
 "NC" 1.003769376117994E-7,
 "GN" 9.99498017070683E-8,
 "CN" 9.781412218341299E-8,
 "AR" 8.60678848033088E-8,
 "TY" 7.667089489922546E-8,
 "CY" 7.154526404245272E-8,
 "GY" 6.684676909041106E-8,
 "RA" 6.300254594783151E-8,
 "RG" 5.232414832955498E-8,
 "RC" 4.890706109170649E-8,
 "YG" 4.847992518697543E-8,
 "AY" 4.847992518697543E-8,
 "CR" 4.741208542514778E-8,
 "CS" 4.036434299708527E-8,
 "MC" 3.929650323525762E-8,
 "CM" 3.780152756869891E-8,
 "YC" 3.780152756869891E-8,
 "SC" 3.673368780687125E-8,
 "AK" 3.652011985450572E-8,
 "RT" 3.630655190214019E-8,
 "WA" 3.459800828321595E-8,
 "TW" 3.331660056902276E-8,
 "AM" 3.310303261665723E-8,
 "GS" 3.310303261665723E-8,
 "GR" 3.267589671192617E-8,
 "TR" 3.203519285482958E-8,
 "KG" 3.182162490246405E-8,
 "YA" 3.054021718827087E-8,
 "SG" 3.054021718827087E-8,
 "MA" 2.989951333117428E-8,
 "AS" 2.968594537880874E-8,
 "GK" 2.883167356934662E-8,
 "KA" 2.840453766461556E-8,
 "TK" 2.712312995042238E-8,
 "KC" 2.690956199805685E-8,
 "ST" 2.669599404569132E-8,
 "KT" 2.58417222362292E-8,
 "MT" 2.349247476020836E-8,
 "CK" 2.30653388554773E-8,
 "WT" 2.285177090311177E-8,
 "SA" 2.263820295074624E-8,
 "GM" 2.178393114128411E-8,
 "AW" 2.157036318891858E-8,
 "CW" 2.114322728418752E-8,
 "MG" 2.114322728418752E-8,
 "TM" 1.986181956999434E-8,
 "WG" 1.793970799870457E-8,
 "TS" 1.623116437978032E-8,
 "WC" 1.580402847504926E-8,
 "GW" 1.452262076085608E-8,
 "RR" 5.339198809138263E-9,
 "MM" 3.84422314257955E-9,
 "SS" 3.630655190214019E-9,
 "KK" 2.989951333117428E-9,
 "YY" 2.989951333117428E-9,
 "YW" 2.989951333117428E-9,
 "KR" 2.776383380751897E-9,
 "SM" 2.562815428386366E-9,
 "WW" 2.349247476020836E-9,
 "SY" 2.349247476020836E-9,
 "MY" 2.349247476020836E-9,
 "WR" 2.349247476020836E-9,
 "RK" 2.135679523655305E-9,
 "WS" 2.135679523655305E-9,
 "WY" 1.922111571289775E-9,
 "RY" 1.922111571289775E-9,
 "KY" 1.922111571289775E-9,
 "YK" 1.922111571289775E-9,
 "YR" 1.922111571289775E-9,
 "SR" 1.922111571289775E-9,
 "KM" 1.708543618924244E-9,
 "SW" 1.708543618924244E-9,
 "KS" 1.708543618924244E-9,
 "MW" 1.708543618924244E-9,
 "RM" 1.708543618924244E-9,
 "BG" 1.494975666558714E-9,
 "CH" 1.494975666558714E-9,
 "MR" 1.494975666558714E-9,
 "RW" 1.494975666558714E-9,
 "GB" 1.494975666558714E-9,
 "RS" 1.281407714193183E-9,
 "VA" 1.281407714193183E-9,
 "BT" 1.281407714193183E-9,
 "CV" 1.281407714193183E-9,
 "WM" 1.281407714193183E-9,
 "SK" 1.281407714193183E-9,
 "YS" 1.281407714193183E-9,
 "CB" 1.281407714193183E-9,
 "AD" 1.067839761827653E-9,
 "MS" 1.067839761827653E-9,
 "KW" 1.067839761827653E-9,
 "VC" 1.067839761827653E-9,
 "AV" 1.067839761827653E-9,
 "HC" 1.067839761827653E-9,
 "CD" 8.542718094621221E-10,
 "TV" 8.542718094621221E-10,
 "DG" 8.542718094621221E-10,
 "HT" 8.542718094621221E-10,
 "TB" 8.542718094621221E-10,
 "GV" 8.542718094621221E-10,
 "VG" 8.542718094621221E-10,
 "YM" 8.542718094621221E-10,
 "BA" 8.542718094621221E-10,
 "DC" 8.542718094621221E-10,
 "HH" 6.407038570965916E-10,
 "GH" 6.407038570965916E-10,
 "WK" 6.407038570965916E-10,
 "SN" 6.407038570965916E-10,
 "VT" 6.407038570965916E-10,
 "AB" 4.271359047310611E-10,
 "DT" 4.271359047310611E-10,
 "TD" 4.271359047310611E-10,
 "YN" 4.271359047310611E-10,
 "DA" 4.271359047310611E-10,
 "GD" 4.271359047310611E-10,
 "NK" 4.271359047310611E-10,
 "MK" 4.271359047310611E-10,
 "HG" 4.271359047310611E-10,
 "NM" 4.271359047310611E-10,
 "BC" 2.135679523655305E-10,
 "MN" 2.135679523655305E-10,
 "VW" 2.135679523655305E-10,
 "KN" 2.135679523655305E-10,
 "VY" 2.135679523655305E-10,
 "NR" 2.135679523655305E-10,
 "AH" 2.135679523655305E-10,
 "DM" 2.135679523655305E-10,
 "MV" 2.135679523655305E-10,
 "NW" 2.135679523655305E-10,
 "HR" 2.135679523655305E-10,
 "NY" 2.135679523655305E-10,
 "YH" 2.135679523655305E-10,
 "BR" 2.135679523655305E-10,
 "TH" 2.135679523655305E-10,
 "WN" 2.135679523655305E-10,
 "HA" 2.135679523655305E-10,
 "RN" 2.135679523655305E-10,
 "VS" 2.135679523655305E-10,
 "WV" 2.135679523655305E-10}


With symmetry

{"GC" 0.1650928245554691,
 "GA" 0.1150176810985652,
 "TC" 0.11501714077164571,
 "CA" 0.11302103730185117,
 "TG" 0.11297672237887124,
 "AT" 0.10850038849932647,
 "TT" 0.0696883223046914,
 "AA" 0.06967938106879766,
 "CC" 0.06552183985818097,
 "GG" 0.06547925910697425,
 "NN" 2.530566667579171E-6,
 "NA" 2.913066870265836E-7,
 "TN" 2.755026585515344E-7,
 "YT" 2.2467348588853815E-7,
 "NG" 2.0246241884252298E-7,
 "NC" 1.981910597952124E-7,
 "AR" 1.4907043075114032E-7,
 "GY" 1.153266942773865E-7,
 "CY" 1.0934679161115162E-7,
 "RC" 9.631914651685427E-8,
 "RG" 8.500004504148115E-8,
 "AY" 7.90201423752463E-8,
 "CS" 7.709803080395653E-8,
 "MC" 7.709803080395653E-8,
 "RT" 6.834174475696978E-8,
 "AK" 6.492465751912127E-8,
 "GS" 6.364324980492809E-8,
 "AM" 6.300254594783151E-8,
 "KG" 6.065329847181067E-8,
 "TW" 5.6168371472134536E-8,
 "WA" 5.616837147213452E-8,
 "TK" 5.2964852186651575E-8,
 "AS" 5.232414832955498E-8,
 "KC" 4.9974900853534154E-8,
 "MT" 4.33542943302027E-8,
 "ST" 4.292715842547164E-8,
 "GM" 4.2927158425471626E-8,
 "CW" 3.694725575923678E-8,
 "WG" 3.246232875956065E-8,
 "RR" 5.339198809138263E-9,
 "YW" 4.912062904407203E-9,
 "KR" 4.9120629044072015E-9,
 "MM" 3.84422314257955E-9,
 "RY" 3.84422314257955E-9,
 "KY" 3.84422314257955E-9,
 "WR" 3.84422314257955E-9,
 "WS" 3.844223142579549E-9,
 "SS" 3.630655190214019E-9,
 "SY" 3.6306551902140187E-9,
 "SM" 3.6306551902140187E-9,
 "MY" 3.203519285482958E-9,
 "RM" 3.203519285482958E-9,
 "SR" 3.203519285482958E-9,
 "KK" 2.989951333117428E-9,
 "YY" 2.989951333117428E-9,
 "BG" 2.989951333117428E-9,
 "KS" 2.989951333117427E-9,
 "MW" 2.989951333117427E-9,
 "CH" 2.5628154283863672E-9,
 "WW" 2.349247476020836E-9,
 "VA" 2.349247476020836E-9,
 "CV" 2.349247476020836E-9,
 "KM" 2.135679523655305E-9,
 "BT" 2.135679523655305E-9,
 "KW" 1.7085436189242446E-9,
 "CD" 1.7085436189242441E-9,
 "GV" 1.7085436189242441E-9,
 "AD" 1.494975666558714E-9,
 "TV" 1.4949756665587137E-9,
 "CB" 1.4949756665587135E-9,
 "DG" 1.2814077141931832E-9,
 "BA" 1.2814077141931832E-9,
 "GH" 1.0678397618276527E-9,
 "HT" 1.0678397618276525E-9,
 "DT" 8.542718094621222E-10,
 "HH" 6.407038570965916E-10,
 "YN" 6.407038570965916E-10,
 "SN" 6.407038570965916E-10,
 "NK" 6.407038570965916E-10,
 "NM" 6.407038570965916E-10,
 "VW" 4.27135904731061E-10,
 "NR" 4.27135904731061E-10,
 "AH" 4.27135904731061E-10,
 "NW" 4.27135904731061E-10,
 "VY" 2.135679523655305E-10,
 "DM" 2.135679523655305E-10,
 "MV" 2.135679523655305E-10,
 "HR" 2.135679523655305E-10,
 "YH" 2.135679523655305E-10,
 "BR" 2.135679523655305E-10,
 "VS" 2.135679523655305E-10}

