;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                              P I P E L I N E                             ;;
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

(ns edu.bc.bio.gaisr.pipeline

  "Primary functions defining the overall pipeline of GAISR."

  (:require [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure-csv.core :as csv]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])

  (:use [clojure.contrib.math :as math]
        [clojure.pprint
         :only [cl-format]]

        [edu.bc.log4clj :only [create-loggers log>]]
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.tools

        [edu.bc.bio.job-config :only [parse-config-file]]

        [edu.bc.bio.gaisr.operon-ctx
         :only [get-region]]
        [edu.bc.bio.gaisr.db-actions
         :only [names->tax]]
        [edu.bc.bio.gaisr.post-db-csv
         :only [+cmsearch-csv-header+]]
        ))




(def blastdb-dir (str "/tmp/blast" (name (gen-kwuid)) "/"))

(defn hit-fields [hit-line]
  (vec (str/split #"," hit-line)))

(defn get-seq-entry [entry-line]
  (let [fields (if (string? entry-line) (hit-fields entry-line) entry-line)
        entry (re-find #"NC_[0-9A-Z]+" (fields 4))
        start (Integer. (fields 5))
        end   (Integer. (fields 6))
        start0 start
        end0   end
        tmp start
        strand (if (< end start) -1 1)
        start (if (< end start) end start)
        end (if (= start end) tmp end)
        range (str start "-" end)]
    [entry (str start0 "-" end0) strand (str entry " " range "/" strand)]))


(defn basic-entries [hits & {form :form :or {form :name-loc}}]
  (keep #(let [[nm loc strand entry] (get-seq-entry %)]
           (when nm
             (case form
                   :name-loc (str nm " " loc)
                   :full entry
                   :full-vec [nm loc strand]
                   :name nm
                   :name-loc-vec [nm loc])))
        hits))


;;; "/home/jsa/Bio/Blasting/b-subtilus-Ls.fna.blast"
(defn hitfile->basic-entries [hitfile & {form :form :or {form :name-loc}}]
  (-> hitfile
      slurp
      ((fn[x](str/split #"\n" x)))
      (basic-entries :form form)))




;;; ----------------------------------------------------------------------
;;; UTR filtering, translation to fasta file, generation of
;;; name/loc/seq tuples


(defn utrs-xform [nm-loc-stg upstream downstream]
  (let [[name loc] (str/split #" " nm-loc-stg)
        [s e] (map #(Integer. %) (str/split #"-" loc))
        strand (if (< e s) -1 1)
        [ns ne] (get-region name strand s upstream downstream)
        tmp ns
        ns (if (< ne ns) ne ns)
        ne (if (= ns ne) tmp ne)
        len (math/abs (- ne ns))]
    [name [ns ne] strand len]))

(defn utrs-filter [entries
                   & {:keys [par upstream downstream]
                      :or {par 10 upstream 500 downstream 25}}]
  (println :--> " " (first entries))
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
             (doall (map #(utrs-xform % upstream downstream) entset)))
           entsets))))


(defn name-seqs [pairs entries]
  (map (fn [nmloc [_ sq]]
         (let [[entry loc] (str/split #" " nmloc)]
           [entry loc sq]))
       entries
       pairs))

(defn utrs-2-entfile-fmt [entries]
  (map (fn[[nm [s e] strand _]] (str nm " " s "-" e "/" strand)) entries))

(defn nlsq-tuples-from-utrs [entries hitfna]
  (let [subdir (let [x (str blastdb-dir (name (gen-kwuid)))]
                 (fs/mkdirs x) x)
        entries (utrs-2-entfile-fmt entries)
        efile (gen-entry-file entries (fs/tempfile "hit-seqs-" ".ent" subdir))
        nlsq-tuples (-> efile
                        entry-file->fasta-file
                        ((fn[x](fs/copy x hitfna) x))
                        slurp
                        ((fn[x](str/split #"\n" x)))
                        ((fn[x] (partition 2 x)))
                        (name-seqs entries))]
    (fs/rm-rf subdir)
    nlsq-tuples))




;;; ----------------------------------------------------------------------
;;; Phylogeny clustering: group and coalesce hits by taxon and cutoff size

(defn entries-at-level [taxon-tuple-map level]
  (group-by second
            (map (fn[m]
                   [(m :name)
                    (let [as (drop 1 (str/split #", " (m :ancestors)))]
                      (if (< level (count as))
                        (nth (reverse as) level)
                        (first as)))])
                 (names->tax taxon-tuple-map))))

(defn phylo-cluster-tuples [taxon-tuple-map level]
  (let [entries (seq (entries-at-level taxon-tuple-map level))]
    (reduce
     (fn[m e]
       (let [mi (reduce (fn [m [n tx]]
                          (assoc m tx (set/union (get m tx #{})
                                                 (set (taxon-tuple-map n)))))
                        {} (val e))]
         (merge-with #(set/union %1 %2) m mi)))
     {} entries)))


(defn sift-clusters [cluster-map cutoff]
  (map #(do [(first %1) (count (second %1)) (second %1)])
       (filter #(> (count (second %1)) cutoff) cluster-map)))

(defn cluster-set-entries [cluster-set]
  (reduce (fn[s e] (conj s (e 0)))
          #{} (apply concat (map #(%1 2) cluster-set))))

(defn remove-tuple-map-entries [taxon-tuple-map entry-set]
  (apply dissoc taxon-tuple-map entry-set))


(defn coalesce-subclusters [clusters]
  (let [tx-cluster-map (group-by first (sort-by second clusters))]
    (map #(do [(key %1)
               (reduce (fn[v [n c sqs :as vx]]
                         (if (empty? v)
                           (conj v vx)
                           (let [[v1n v1c v1sqs :as v1] (first v)]
                             ;;(prn (str "***" c " " v1c) (map second v))
                             (if (and (>= (+ v1c c) 200) (>= c 50))
                               (conj v vx)
                               (let [comb-sqs (set/union sqs v1sqs)
                                     comb-cnt (count comb-sqs)
                                     newv [n comb-cnt comb-sqs]]
                                 (conj (drop 1 v) newv))))))
                       ()
                       (reverse (val %1)))]) tx-cluster-map)))

(defn flatten-to-subclusters [clusters]
  (apply concat (map (fn[[k v]] v) clusters)))

;;; Debugging/Checking vars.
(defparameter *tmap* (atom {}))
(defparameter *tclusts* (atom ()))

(defn group-by-taxon [taxon-tuple-map & {cutoff :cutoff :or {cutoff 100}}]
  (loop [clusters ()
         todo-map taxon-tuple-map
         level 1
         cutoff cutoff]
    ;;(swap! *tmap* (fn[_] todo-map))
    ;;(swap! *tclusts* (fn[_] clusters))
    (cond
     (or (empty? todo-map) (< cutoff 3))
     (flatten-to-subclusters (coalesce-subclusters clusters))

     (> level 20)
     (recur clusters todo-map 1 (math/floor (/ cutoff 2)))

     :else
      (let [cluster-map (phylo-cluster-tuples todo-map level)
            cutoff-clusters (sift-clusters cluster-map cutoff)
            entry-set (cluster-set-entries cutoff-clusters)]
        (recur (set/union clusters cutoff-clusters)
               (remove-tuple-map-entries todo-map entry-set)
               (inc level)
               cutoff)))))


(defn get-taxon-cluster-map [hitfile]
  (let [nlsq-tuples
        (-> hitfile
            hitfile->basic-entries
            utrs-filter
            (nlsq-tuples-from-utrs (fs/replace-type hitfile ".hitfna")))]
    (reduce (fn[m [nm rng sq]]
              (assoc m nm (conj (get m nm []) [nm rng sq])))
            {} nlsq-tuples)))


(defn phylo-clusters [hitfile]
  (-> hitfile get-taxon-cluster-map group-by-taxon))




;;; ----------------------------------------------------------------------
;;;
;;; Redundancy computation/filtering Full Levenshtein global or ngram
;;; closeness or statistical/levenshtein mix via cdhit.


(defn edist [[nq lq qseq] [nms ls sseq]]
  (let [qlen (count qseq)
        slen (count sseq)
        %70-qlen (* 0.7 qlen)]
    (if (<= slen %70-qlen)
      [:keep [slen %70-qlen] 0.0 0 [nms ls sseq]]
      (let [dist (levenshtein qseq sseq)
            %ident (float (/ (math/abs (- dist qlen)) qlen))]
        [(if (>= %ident 0.9) :toss :keep)
         [slen %70-qlen] %ident dist [nms ls sseq]]))))

(defn ngram-closeness [[nq lq qseq] [nms ls sseq]
                       & {c :c n :n :or {c 0.85 n 4}}]
  (let [qlen (count qseq)
        slen (count sseq)
        %70-qlen (* 0.7 qlen)]
    (if (<= slen %70-qlen)
      [:keep [slen %70-qlen] 0.0 0 [nms ls sseq]]
      (let [%ident (float (ngram-compare qseq sseq :n n :uc? false))]
        [(if (>= %ident c) :toss :keep)
         [slen %70-qlen] %ident 0 [nms ls sseq]]))))


(defn seq-same? [qseq test-set & {cmpfn :cmpfn :or {cmpfn edist}}]
  (some #(= :toss (first (cmpfn qseq %)))
        test-set))

(defn pseq-same? [qseq test-set & {q :q cmpfn :cmpfn :or {q 10 cmpfn edist}}]
  (loop [curset test-set
         chunk (take q curset)]
    (cond
     (empty? chunk) false
     (some #{:toss} (pmap #(first (cmpfn qseq %)) chunk)) true
     :else (let [nextset (drop q curset)]
             (recur nextset (take q nextset))))))


(defn get-keeper-set [nlsq-tuples & {cmpfn :cmpfn :or {cmpfn edist}}]
  (let [tuples (sort-by #(% 2) (fn[l r] (> (count l) (count r))) nlsq-tuples)]
    (binding [seq-same? (if (< (count tuples) 100) seq-same? pseq-same?)]
      (loop [keepers []
             todo tuples]
        (if (empty? todo)
          keepers
          (let [rep (first todo)]
            (recur (if (not (seq-same? rep keepers :cmpfn cmpfn))
                     (conj keepers rep)
                     keepers)
                   (rest todo))))))))


(defparameter *clusnr* "ClusNR")

(defn get-candidates [hit-file & {:keys [cmpfn dir]
                                  :or {cmpfn cd-hit-est dir nil}}]
  (let [hit-file (fs/fullpath hit-file)
        dir (fs/fullpath (fs/join (if dir dir (fs/dirname hit-file)) *clusnr*))
        dir? (or (fs/exists? dir) (fs/mkdir dir))
        s2- #(str/replace-re #" " "-" %)
        tuple-clusters (phylo-clusters hit-file)]
    (when (not dir?)
      (raise :type :dir-notexist :path dir))
    (if (= cmpfn cd-hit-est)
      (doseq [[tx cnt nrsq-tuples] tuple-clusters]
        (cd-hit-est nrsq-tuples (fs/join dir (str (s2- tx) "-" cnt ".fna"))))
      (let [clusters (map (fn[[tx cnt nrsq-tuples]]
                            [tx cnt (get-keeper-set nrsq-tuples :cmpfn cmpfn)])
                          tuple-clusters)]
        (doseq [[tx cnt nrsq-tuples] clusters]
          (nms-sqs->fasta-file
           (map (fn[[nm rng sq]] [(str nm ":" rng) sq]) nrsq-tuples)
           (fs/join dir (str (s2- tx) "-" cnt ".fna"))))
        clusters))))




;;; ----------------------------------------------------------------------
;;;

(defn summary-map [motif-sto-file]
  (let [m (into {}
                (map #(let [[x y] (str/split #"=" (str/trim %))]
                        [(keyword (str/replace-re #"_" "-" (str/lower-case x)))
                         (Float. y)])
                     (str/split
                      #"\t" (runx "summarize" "-w" motif-sto-file))))
        species-cnt (count
                     (frequencies
                      (map #(second (re-find #"#=GS\s+(\S+):" %))
                           (filter #(re-find #"DE" %)
                                   (str/split #"\n" (slurp motif-sto-file))))))]
    (assoc m :species-cnt species-cnt :motif motif-sto-file)))

(defn rank-em [summary-maps]
  (sort-by
   #(% :rank) >
   (map
    (fn[sm]
      (let [sid (sm :seq-id) ; % seq identical
            sid (if (<= sid 0) 1 sid)]
        (assoc sm
          :rank
          (* (sm :species-cnt) (math/sqrt (* (+ (sm :conserved-pos) 0.2)
                                         (/ (sm :bp) sid)))
             (inc (java.lang.Math/log (/ (sm :num) (sm :species-cnt))))))))
    summary-maps)))

;;;ord <- order(test.data[,"Rank.index"], decreasing=T);


(defn cmf-post-process [motif-sto-filespec
                        & {lt :lt ut :ut w :w
                           :or {lt 10 ut nil s true w nil}}]
  (let [infile (fs/fullpath motif-sto-filespec)
        outfile (str infile ".filtered")
        cmfpath (get-tool-path :cmfinder)
        cmfiltercmd (str cmfpath "filter.pl")
        cmdargs ["-s" "-lt" (str lt)]
        cmdargs (if ut (conj cmdargs "-ut" (str ut)) cmdargs)
        cmdargs (conj cmdargs infile outfile)]
    (assert-tools-exist [cmfiltercmd])
    (runx cmfiltercmd cmdargs)
    outfile))


(defn cmbuild [motif-stofile]
  (let [infernal-path (get-tool-path :infernal)
        cmbuildcmd (str infernal-path "cmbuild")
        motif-stofile (fs/fullpath motif-stofile)
        cmfile (-> motif-stofile
                   (#(str/split #"\." %))
                   (#(str (first %) "."
                          (last (if (> (count %) 2) (butlast %) %))))
                   (str ".cm"))]
    (assert-tools-exist [cmbuildcmd])
    (runx cmbuildcmd "-F" cmfile motif-stofile)
    cmfile))


(defn cmcalibrate [cmfile & {par :par  :or {par 3}}]
  (let [infernal-path (get-tool-path :infernal)
        cmcalibratecmd (str infernal-path "cmcalibrate")
        cmfile (fs/fullpath cmfile)
        mpirun "mpirun"
        cmdargs ["-np" (str par)
                 cmcalibratecmd "--mpi" cmfile]]
    (assert-tools-exist [cmcalibratecmd])
    (runx mpirun cmdargs)
    cmfile))


(defn cmsearch [cmfile fna-seqalign-file outfile
                & {par :par eval :eval :or {par 3 eval 1.0}}]
  (let [infernal-path (get-tool-path :infernal)
        cmsearchcmd (str infernal-path "cmsearch")
        cmfile (fs/fullpath cmfile)
        alignfile (fs/fullpath fna-seqalign-file)
        outfile (fs/fullpath outfile)
        mpirun "mpirun"
        cmdargs ["-np" (str par)
                 cmsearchcmd "--mpi" "-E" (str eval) cmfile alignfile]]
    (assert-tools-exist [cmsearchcmd])
    (io/with-out-writer outfile
      (print (runx mpirun cmdargs)))
    outfile))




;;; ----------------------------------------------------------------------
;;;
;;; CMSearch output filter and CSV representation.  Parses cmsearch
;;; output files, for each plus/minus hit take the better of the two
;;; (in the often case of only one, just take that).


(defn get-generic-sto-locs [stofile-lazyseq]
  (reduce (fn[m s-loc]
            (let [[nm loc] (str/split #"/" s-loc)
                  [nm v] (str/split #"\." nm)
                  [s e] (map #(Integer. %) (str/split #"-" loc))
                  x s
                  s (if (< s e) s e)
                  e (if (= s e) x e)]
              (assoc m nm (conj (get m nm []) [s e]))))
          {} (map #(first (str/split #"\s+" %))
                  (filter #(re-find #"^N(C|S|Z)" %) stofile-lazyseq))))

(defn get-cmfinder-sto-locs [stofile-lazyseq]
  (reduce (fn [m s-loc]
            (let [[nm loc] (str/split #":" s-loc )
                  [loc strand] (str/split #"/" loc)
                  [s e] (str/split #"-" loc)
                  s (Integer. s)
                  e (Integer. e)]
              (assoc m nm (conj (get m nm []) [s e]))))
          {} (map #(second (str/split #"\s+" %))
                  (filter #(re-find #"DE" %) stofile-lazyseq))))

(defn get-infernal-sto-locs [stofile-lazyseq]
  (reduce (fn [m s-loc]
            (let [[nm loc] (str/split #"(/c|/)" s-loc)
                  [nm v] (str/split #"\." nm)
                  [s e] (str/split #"-" loc)
                  s (Integer. s)
                  e (Integer. e)
                  x s
                  s (if (< s e) s e)
                  e (if (= s e) x e)]
              (assoc m nm (conj (get m nm []) [s e]))))
          {} (map #(first (str/split #"\s+" %))
               (filter #(re-find #"^N(C|S|Z)" %) stofile-lazyseq))))

(defn get-sto-seq-locs [stofile]
  (let [stofile-lazyseq (io/read-lines stofile)
        fmt-line (second stofile-lazyseq)
        file-fmt (cond
                  (re-find #"CMfinder" fmt-line) :cmfinder
                  (re-find #"Infernal" fmt-line) :infernal
                  (re-find #"ORIGINAL_MOTIF" fmt-line) :generic
                  :else (raise :type :unknown-sto-fmt
                               :args [stofile fmt-line]))
        rem-file (drop-until #(not (.startsWith % "#")) stofile-lazyseq)]
    (case file-fmt
          :cmfinder (get-cmfinder-sto-locs rem-file)
          :infernal (get-infernal-sto-locs rem-file)
          :generic  (get-generic-sto-locs rem-file))))


(defn get-cm-stofile [cmfile]
  (last (str/split #" " (first (drop-until
                                #(re-find #"^BCOM" %)
                                (str/split #"\n" (slurp cmfile)))))))

(defn ent-name [ent-line]
  (if (= (str/take 3 ent-line) ">gi")
    (first (re-find #"N(C|S|Z)_[0-9A-Z]+" ent-line))
    (first (entry-parts ent-line))))

(defn ent-loc [ent-line]
  (if (= (str/take 3 ent-line) ">gi")
    (let [l (first (re-find #":(.|)[0-9]+-[0-9]+" ent-line))]
      (when l (subs l (first (pos-any "0123456789" l)))))
    (->> ent-line
         entry-parts
         ((fn[[nm [s e] sd]]
            (let [[s e] (if (= sd "-1") [e s] [s e])]
              (str s "-" e)))))))

(defn build-hitseq-map [hitfile]
  (reduce (fn[m [gi sq]]
            (let [nc (ent-name gi)
                  k (if nc (str nc ":" (ent-loc gi)) gi)]
              (assoc m k sq)))
          {} (partition 2 (io/read-lines (io/file-str hitfile)))))


(defn cmsearch-group-hits [cmsearch-out]
  (let [file-content (str/split #"\n" (slurp cmsearch-out))
        cmdline (first (drop-until #(re-find #"^# command" %) file-content))
        cmfile (subs (re-find #" /[A-Za-z0-9_\.\-/]+\.cm" cmdline) 1)
        hitfile (subs (first (re-find #" /[A-Za-z0-9_\.\-/]+\.(hitfna|fa|fna)"
                                      cmdline)) 1)
        stofile (get-cm-stofile cmfile)
        lines (drop-until #(re-find #"^>" %) file-content)
        parts (partition-by #(if (or (= % "#") (re-find #"^>" %)) :x :y)
                            lines)]
    [(get-sto-seq-locs stofile)
     (reduce (fn[m [k v]]
               (if (= (first k) "#")
                 m
                 (let [k (first k)
                       nc (ent-name k)
                       k (if nc (str nc ":" (ent-loc k)) k)
                       v (keep #(when (not= "" %) (str/trim %)) v)]
                   (assoc m k v))))
             {} (partition 2 parts))
     (build-hitseq-map hitfile)
     stofile]))


(defn remove-pre&sufix-locs [seq-stg]
  (str/replace-re #"(^[0-9]+ | [0-9]+$)" "" seq-stg))

(defn cmsearch-hit-info [hit-val]
  (map (fn[[hit-strand query-tgt sc-ev-p-gc & tail]]
         (flatten [(if (= "Plus" (subs hit-strand 0 4)) 1 -1)
                   (let [[s _ e] (drop 7 (str/split #" " query-tgt))]
                     [(Integer. s) (Integer. e)])
                   (map #(second (str/split #" = *" %))
                        (str/split #", " sc-ev-p-gc))
                   (let [x (partition 4 tail)
                         x (apply map (fn[& tail]
                                        (apply str (map remove-pre&sufix-locs
                                                        tail)))
                                  x)
                         struct (first x)
                         sq (str/replace-re #"(^[0-9]+ | [0-9]+$)" "" (last x))]
                          [struct sq])]))
       (map flatten
            (partition 2 (partition-by #(if (re-find #"strand" %) :x :y)
                                       hit-val)))))

(defn get-hit-loc [seq-start seq-end hit-start hit-end]
  (let [strand (if (> seq-start seq-end) -1 1)]
    [(+ seq-start (* strand (dec hit-start)))
     (+ seq-start (* strand hit-end))]))

(defn get-orig-seq-full-hit
  "Return full original sequence of the hit location given by hit-start S and
   hit-end E over original sequence ORIG-SEQ.  If strand is minus return the
   reverse compliment."
  [orig-seq s e]
  (let [[s e st] (if (> s e) [e s -1] [s e 1])
        twiddle (if (= st -1) reverse-compliment identity)]
    (str/replace-re #"T" "U" (twiddle (subs orig-seq (dec s) e)))))

(defn cmsearch-hit-parts [[h v] hit-seq-map]
  (let [[nm loc] (str/split #":" h)
        orig-seq (hit-seq-map h)
        [s e] (if loc (vec (str/split #"-" loc)) ["1" (str (count orig-seq))])
        s (Integer. s)
        e (Integer. e)
        info (reduce
              (fn[cur nxt] ; Keep the one with best Evalue
                (if (< (Float. (nth nxt 4)) (Float. (nth cur 4)))
                  nxt
                  cur))
              [0 1 2 3 100000.0] (cmsearch-hit-info v))
        [_ hit-start hit-end & tail] info
        [ns ne] (get-hit-loc s e hit-start hit-end)]
    (vec (flatten [nm s e ns ne info
                   (get-orig-seq-full-hit orig-seq hit-start hit-end)]))))


(defn group-hit-parts [sto-loc-map hit-parts]
  (group-by
   (fn[[nm s e & tail]]
     ;;(prn (= @*dbg-hit-parts* hit-parts) nm (count (sto-loc-map nm)) s e)
     (let [locs (sort (sto-loc-map nm))
           [s e] (sort [s e])]
       (if (or (empty? locs)
               (empty? (keep (fn[[ls le :as v]]
                               (when (not (or (< ls (+ s 10) le)
                                              (< ls (- e 10) le)))
                                 v))
                             locs)))
         :good :bad)))
   hit-parts))

(defn remove-dups
  "Hits are relative in cmsearch out (sto file).  The original hit
   start/end are relative to the relative sequence in the hit
   information and the hit-loc-map we create uses these coordinates to
   attempt duplication removal, but these relative coordinates when
   changed to absolute coordinates for the entire genome seq can be
   identical.  We don't have absolute coordinates when building the
   original map, but we do have them after CMSEARCH-HIT-PARTS.  We
   maybe could remove these duplicates earlier (and not make yet
   another pass!), but this at least has the advantage of being
   transparent.

   hit-parts is a seq of vecs [nm os oe abs abe & tail], where os and
   oe are the original relative hit start/end and abs and abe are the
   computed absolute hit start/end.

   spread is positive integer or 0 which defines the max distance
   between the abs and abe of any two hit-parts for them to be
   considered the same hit; defaults to 10 bases.

   stofile is the originating sto file for the cmsearch output that is
   now encoded in hit-parts.  We ensure there are no entries in the
   output whose coordinates are already within spread of an existing
   entry in the originating sto.

   Filter hit-parts so only one [nm abs abe] identified (keyed) entry
   remains.
  "
  [hit-parts spread stofile]

  (let [sto-map (reduce (fn[M e]
                          (let [[nm coord] (->> e entry-parts (take 2))]
                            (assoc M nm coord)))
                        {} (read-seqs stofile :info :name))
        in-spread? (fn [x y] (<= (abs (- x y)) spread))
        dup? (fn[m [nm abst abend :as k]]
               (let [[s e :as coord] (m k)
                     [stost stoend :as sto-coord] (sto-map nm)]
                 (or (and coord
                          (in-spread? s abst)
                          (in-spread? e abend))
                     (and sto-coord
                          (let [[abst abend]
                                (if (< abend abst) [abend abst] [abst abend])]
                            (in-spread? stost abst)
                            (in-spread? stoend abend))))))]

    (second ; just return filtered seq, not the map used...
     (reduce (fn [[m s] v]
               (let [[nm os oe abst abend & tail] v
                     k [nm abst abend]]
                 (if (dup? m k)
                   [m s]
                   [(assoc m k [abst abend]) (conj s v)])))
             [{} []] hit-parts))))


(defn cmsearch-out-csv [cmsearch-out & {:keys [spread] :or {spread 50}}]
  (if (fs/empty? cmsearch-out)
    [[] []]
    (let [csv-file (str/replace-re #"\.out$" ".csv" cmsearch-out)
          dup-file (str/replace-re #"\.out$" ".dup.csv" cmsearch-out)

          [sto-loc-map hit-map hit-seq-map stofile]
          (cmsearch-group-hits cmsearch-out)

          hit-parts (map #(cmsearch-hit-parts % hit-seq-map) hit-map)
          groups (group-hit-parts sto-loc-map hit-parts)

          good (map #(csv/csv-to-stg (map str %))
                    (remove-dups (:good groups) spread stofile))
          dups (map #(csv/csv-to-stg (map str %))
                    (:bad groups))]

      (io/with-out-writer (io/output-stream csv-file)
        (println (str +cmsearch-csv-header+ "," stofile))
        (doseq [x good] (println x)))
      (io/with-out-writer (io/output-stream dup-file)
        (println +cmsearch-csv-header+)
        (doseq [x dups] (println x)))
      [good dups])))


(defn gen-cmsearch-csvs [cmsearch-out-dir & {:keys [spread] :or {spread 50}}]
  (let [base cmsearch-out-dir]
    (doseq [x (filter #(re-find #"\.cmsearch\.out$" %)
                      (sort (map #(fs/join base %) (fs/listdir base))))]
      (cmsearch-out-csv x :spread spread))))




;;; ----------------------------------------------------------------------
;;;
;;; Upper layer pipeline threading and running.

(defn create-pipeline-working-area [ClusNR-dir]
  (let [base (fs/dirname ClusNR-dir)
        clus-files (sort (filter #(not (re-find #"clstr" %))
                                 (fs/listdir ClusNR-dir)))
        clus-dirnms&fs (map #(do [(first (str/split #"\." %)) %])
                            clus-files)]
    (doseq [[dnm fs] clus-dirnms&fs]
      (let [ndir (fs/join base dnm)]
        (when (not (fs/exists? ndir))
          (fs/mkdir ndir))
        (fs/copy (fs/join ClusNR-dir fs) (fs/join ndir fs))))
    [base clus-dirnms&fs]))


(defn gen-hit-out-filespec [cm-filespec & {hitfna :hitfna :or {hitfna ""}}]
  (let [hitpart (fs/replace-type (fs/basename hitfna) ".")]
    (fs/join
     (fs/dirname cm-filespec)
     (str/join
      "." (conj (vec (butlast (str/split #"\." (fs/basename cm-filespec))))
                (str hitpart "cmsearch.out"))))))

(defn core-processing [hit-fna base clus-dirnms&fnms]
  (doall
   (map (fn [[dnm fnm]]
          (-> (fs/join base dnm fnm)
              cmfinder*
              ((fn[mstos](map #(cmf-post-process %) mstos)))
              ((fn[fmstos](map #(cmbuild %) fmstos)))
              ((fn[cms](pmap #(cmcalibrate %) cms)))
              ((fn[cms](pmap #(cmsearch % hit-fna (gen-hit-out-filespec %))
                             cms)))))
        clus-dirnms&fnms)))


(defn process-clusters [hit-fna base clus-dirnms&fnms
                        & {par :par :or {par 3}}]
  (let [sz (floor (/ (count clus-dirnms&fnms) par))
        wsets (partition sz sz {} clus-dirnms&fnms)]
    (pmap #(core-processing hit-fna base %) wsets)))


(defn run-pipeline-front [selections
                          & {:keys [ev wordsize] :or {ev 10}}]
  (let [selections-fna (get-selection-fna selections)
        blaster (blastpgm selections-fna)
        wdsz (if wordsize wordsize (if (= blaster tblastn) 4 8))
        hit-file (blaster selections-fna :word-size wdsz :evalue ev)
        hitfna (fs/replace-type hit-file ".hitfna")
        clusters (get-candidates hit-file)]
    [hitfna hit-file]))

(defn run-pipeline-full [selections
                         & {:keys [ev wordsize] :or {ev 10}}]
  (let [selections-fna (get-selection-fna selections)
        blaster (blastpgm selections-fna)
        wdsz (if (= blaster tblastn) 4 8)
        hit-file (blaster selections-fna :word-size wdsz :evalue ev)
        hitfna (fs/replace-type hit-file ".hitfna")
        clusters (get-candidates hit-file)
        clusnr-dir (fs/join (fs/dirname hit-file) *clusnr*)
        [base dir-info] (create-pipeline-working-area clusnr-dir)]
    (doall (process-clusters hitfna base (take 6 dir-info)))))


(defn directory-cms [directory]
  (fs/directory-files directory ".cm"))

(defn directory-hitfnas [directory]
  (fs/directory-files directory ".hitfna"))


(defn cms->calibrated-cms [cms & {par :par :or {par 4}}]
  (loop [cnms (ensure-vec cms)
         results []]
    (let [nextgrp (take 3 cms)]
      (if (empty? nextgrp)
        (flatten results)
        (recur
         (drop 3 cms)
         (conj results
               (doall (pmap #(cmcalibrate % :par par) nextgrp))))))))

(defn mostos->calibrated-cms [mostos & {par :par :or {par 4}}]
  (loop [mostos (ensure-vec mostos)
         results []]
    (let [nextgrp (take 3 mostos)]
      (if (empty? nextgrp)
        (flatten results)
        (recur
         (drop 3 mostos)
         (conj results
               (-> nextgrp
                   ((fn[mstos]
                      (map #(cmbuild %) mstos)))
                   ((fn[cms]
                      (doall (pmap #(cmcalibrate % :par par) cms)))))))))))


;;; (defn mostos&hitfile->cmsearch-out [mostos hifile]
;;;   (let [hit-fna (fs/replace-type hitfile ".hitfna")]
;;;     (-> hitfile
;;;         hitfile->basic-entries
;;;         utrs-filter
;;;         (nlsq-tuples-from-utrs hit-fna))
;;;     ???))


(defn mostos&hitfna->cmsearch-out [mostos hit-fna]
  (-> (ensure-vec mostos)
      ((fn[mstos](map #(cmbuild %) mstos)))
      ((fn[cms](doall (pmap #(cmcalibrate %) cms))))
      ((fn[cms](pmap #(cmsearch % hit-fna (gen-hit-out-filespec %))
                     cms)))))


(defn cms&hitfna->cmsearch-out
  [cms hit-fna & {par :par eval :eval :or {par 3 eval 1.0}}]
  (loop [cms (ensure-vec cms)
         results []]
    (let [nextgrp (take par cms)]
      (if (empty? nextgrp)
        (flatten results)
        (recur
         (drop par cms)
         (conj results
               (doall (pmap
                       #(cmsearch
                         % hit-fna (gen-hit-out-filespec % :hitfna hit-fna)
                         :eval eval)
                       nextgrp))))))))

(defn cms&hitfnas->cmsearch-out
  [cms hit-fnas &
   {cmpar :cmpar par :par eval :eval :or {cmpar 3 par false eval 1.0}}]
  (let [mapper (if par pmap map)]
    (mapper #(cms&hitfna->cmsearch-out cms % :par cmpar :eval eval)
            (ensure-vec hit-fnas))))


(defn cms&hitfile->cmsearch-out
  [cms hitfile & {upstream :upstream downstream :downstream
                  :or {upstream 500 downstream 25}}]
  (let [hit-fna (fs/replace-type hitfile ".hitfna")
        utrs-filter #(utrs-filter % :upstream upstream :downstream downstream)]
    (-> hitfile
        hitfile->basic-entries
        utrs-filter
        (nlsq-tuples-from-utrs hit-fna))
    (pmap #(cmsearch % hit-fna (gen-hit-out-filespec % :hitfna hit-fna))
          (ensure-vec cms))))



(defn mostos&hitfiles->cmsearch-out [mostos hitfiles]
  (let [cms (mostos->calibrated-cms mostos)]
    (pmap #(cms&hitfile->cmsearch-out cms %)
          (let [base "/data2/Bio/ECRibLeaders/FastaFiles"]
            (keep #(when (re-find #"blast$" %) (fs/join base %))
                  (fs/listdir base))))))


;;; Running pipeline via config files and directives

(comment

  (use '[edu.bc.bio.job-config :only (parse-config-file)])
  (fs/directory-files "/home/kaila/Bio/STOfiles/STO-3-100" "sto")
  (fs/exists?
   ((parse-config-file "/home/kaila/Bio/STOfiles/STO-3txt") :gen-csvs))
)


(defn chk-stos [stofiles]
  (let [chk-info (filter (fn [f]
                           (let [r (check-sto f :printem false)]
                             (when (not= r :good) [f r])))
                         stofiles)]
    (when (seq chk-info)
      (raise :type :badsto :chk-info chk-info))))


(defn do-cmbuild-calibrate [stofiles]
  (doall (mostos->calibrated-cms stofiles)))

(defn do-cmbuild [stofiles]
  {:pre [(seq stofiles)]}
  (doall (map #(cmbuild %) stofiles)))

(defn do-calibrate [cms cmdir]
  {:pre [(seq cms)]}
  (let [cmfiles (reduce (fn[v cmfspec]
                          (concat v (fs/re-directory-files cmdir cmfspec)))
                        [] cms)]
    (doall (cms->calibrated-cms cmfiles))))


(defn get-cmsearch-groups [cmdir cmsearchs]
  (map (fn[[hf cmfspecs]]
         [hf (flatten (map #(fs/re-directory-files cmdir %) cmfspecs))])
       cmsearchs)) ; map of hitf->[regex-cm-filespecs]

(defn do-cmsearch [hfs-cmss eval]
  (doall (map (fn[[hf cms]]
                (cms&hitfna->cmsearch-out cms hf :eval eval))
              hfs-cmss)))


(defn run-config-job [job-config-file & {:keys [eval] :or {eval 100.0}}]
  (let [config (parse-config-file job-config-file)
        stodir (config :stodir)
        cmdir (config :cmdir)
        hitdir (config :hitfile-dir)
        stos (or (seq (config :cmbuilds))
                 (seq (fs/directory-files stodir "sto")))
        cms  (or (seq (config :calibrates))
                 (seq (fs/directory-files cmdir "cm")))]

    (when (and (config :check-sto) (config :cmbuild))
      (chk-stos stos))

    (cond
     (and (config :cmbuild) (config :cmcalibrate))
     (do-cmbuild-calibrate stos)

     (config :cmbuild)
     (do-cmbuild stos)

     ;; Note, can't get here unless we have _only_ calibrate directive
     ;; (no build), and so cms can't legitimately be nil.  That is a
     ;; precondition check in do-calibrate
     (config :cmcalibrate)
     (do-calibrate cms cmdir))

    ;; Map over all cmsearch requests
    (let [eval (or (config :eval) eval)
          cmsearchs (config :cmsearchs)
          hfs-cmss (get-cmsearch-groups cmdir cmsearchs)
          cmouts (if (and cmsearchs (seq hfs-cmss))
                   (flatten (do-cmsearch hfs-cmss eval))
                   (fs/directory-files cmdir "cmsearch.out"))]

      ;; Generate csvs for all cmsearch.out's in the cmdir.  If the
      ;; csvdir does not exist, generate it.  Place all generated csvs
      ;; in to csvdir
      (when-let [csvdir (config :gen-csvs)]
        (when (not (fs/exists? csvdir)) (fs/mkdir csvdir))
        (when (seq cmouts)
          (doseq [cmout cmouts] (cmsearch-out-csv cmout))
          (let [csvs (fs/directory-files cmdir "csv")]
            (doseq [csv csvs]
              (let [fname (fs/basename csv)]
                (fs/rename csv (fs/join csvdir fname)))))))
      :good)))

(defn run-config-job-checked
  "Run a configuration with catch and print for any exceptions"
  [job-config-file & {:keys [eval printem] :or {eval 100.0 printem true}}]
  (let [result (catch-all (run-config-job job-config-file :eval eval))]
    (if printem
      (prn result)
      result)))





;;;---------------------- testing - workin progress - messes ------------------


;;; (defn run-pipeline
;;;   [input &
;;;    {form :form steps :steps
;;;     :or {form :full
;;;          steps [hitfile->basic-entries
;;;                 utrs-filter
;;;                 (fn [entries]
;;;                   (nlsq-tuples-from-utrs
;;;                    entries (fs/replace-type input ".hitfna")))]}}]
;;;   (if (= form :full)
;;;     (run-pipeline-full input)
;;;     (let
;;;         )))


;;; (defn pipe-test [hit-fna base [dnm fnm]]
;;;   (-> (fs/join base dnm fnm)
;;;       cmfinder*
;;;       ((fn[mstos](map #(cmf-post-process %) mstos)))
;;;       ((fn[fmstos](map #(cmbuild %) fmstos)))
;;;       ((fn[cms](doall (pmap #(cmcalibrate %) cms))))
;;;       ((fn[cms](pmap #(cmsearch % hit-fna (gen-hit-out-filespec %))
;;;                      cms)))))







;;; (tblastn "/home/jsa/Bio/Blasting/b-subtilus-Ls.fna")
;;; Once through the full pipeline.  Set up a job network embodying a
;;; full pipeline:
;;;


;;; (loop [result (-> selections  ; Selections uploaded or ???
;;;                   blast       ; tblastn or blastn (by selections alphabet)
;;;                   get-candidates ; phylo cluster and redundant filter
;;;                   cmfinder    ; build motif sto files
;;;                   cmfilter    ; filter motif sto files for redundancy
;;;                   cmbuild     ; motif sto files -> infernal cm
;;;                   cmcalibrate ; calibrate cm
;;;                   cmsearch)   ; blast hits + cm --> new alignments
;;;        results #{}]
;;;   (if (no-change result (results :last))
;;;     result ; if change less than delta, return result
;;;     (recur
;;;      (-> result cmbuild cmcalibrate cmsearch)
;;;      (assoc last :last result (gen-kwuid) result))))

;;;filter.pl -s -lt 10 Firmicutes-109.fna.mofif-sto.h1.1 Firmicutes-109-filtered.motif-sto.h1.1
;;;cmbuild Firmicutes-109.cm Firmicutes-109-filtered.motif-sto.h1.1
;;;mpirun -np 4 cmcalibrate --mpi Firmicutes-109.cm
;;;mpirun -np 4 cmsearch --mpi Firmicutes-109.cm ../b-subtilus-Ls.hitfna

;;; (rank-em
;;;  (map summary-map
;;;       (map #(fs/join "/home/jsa/Bio/Pipeline1/Actinobacteridae-103" %)
;;;            (filter #(re-find #"filtered" %)
;;;                    (fs/listdir
;;;                     "/home/jsa/Bio/Pipeline1/Actinobacteridae-103")))))

;;; (map #(do [(% :seq-id) (% :weight)]) *1)


;;; (def *ents* (catch-all (-> "/home/jsa/Bio/Blasting/b-subtilus-Ls.fna.blast"
;;;                            hitfile->basic-entries)))

;;; (def *ents* (time (utrs-filter canned/*ents*)))




;;; (def *x* (catch-all
;;;           (let [hitfile "/data2/Bio/Paper/FastaFiles/Assortprots2.fna.blast"
;;;                 hit-fna (fs/replace-type hitfile ".hitfna")]
;;;             (-> hitfile
;;;                 hitfile->basic-entries
;;;                 utrs-filter
;;;                 (nlsq-tuples-from-utrs hit-fna)))))







