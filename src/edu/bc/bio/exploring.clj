
(in-ns 'edu.bc.bio.gaisr.pipeline)


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


(defn build-hitseq-map [hitfile]
  (reduce (fn[m [gi sq]]
            (let [nc (first (re-find #"N(C|S|Z)_[0-9A-Z]+" gi))
                  k (if nc (str nc ":" (re-find #"[0-9]+-[0-9]+" gi)) gi)]
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
                       nc (first (re-find #"N(C|S|Z)_[0-9A-Z]+" k))
                       k (if nc (str nc ":" (re-find #"[0-9]+-[0-9]+$" k)) k)
                       v (keep #(when (not= "" %) (str/trim %)) v)]
                   (assoc m k v))))
             {} (partition 2 parts))
     (build-hitseq-map hitfile)]))


(catch-all
 (doseq [d (butlast
            (drop 2 (sort (directory-files "/data2/Bio/MetaG1/MoStos" ""))))]
   (gen-cmsearch-csvs d)))