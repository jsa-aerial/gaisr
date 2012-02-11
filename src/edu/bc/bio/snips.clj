
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