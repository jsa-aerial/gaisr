(in-ns 'edu.bc.bio.gaisr.pipeline)


(create-loggers
 [[:logger "metag1-runs"
   {:level :info :appenders "metag1"}]
  [:appender "metag1"
   {:type :rolling
    :filespec "./Logs/metag1-runs.log"
    :max-version 10
    :max-size 1000000}]
  ])

(log> "metag1-runs" :info "hi there - you working?")


(let [base  "/data2/Bio/MetaG1/MoStos"
      dirs (filter  #(fs/directory?
                      (fs/join base %)) (fs/listdir base))]
  (doseq [d dirs]
    (let [f (str d ".txt")
          sto (str d ".sto")]
      (regroup-sto-fasta-file (fs/join base d f) (fs/join base d sto)))))


(def *calibrations*
     (future (let [base  "/data2/Bio/MetaG1/MoStos"
                   dirs (filter  #(fs/directory?
                                   (fs/join base %)) (fs/listdir base))]
               (doseq [d dirs]
                 (let [sto (str d ".sto")]
                   (mostos->calibrated-cms (fs/join base d sto)))))))

(future-done? *calibrations*)
(future-done? *metag1-cmsearch-RF00162*)

(def *metag1-cmsearch-RF00167*
     (future
      (let [*out* edu.bc.bio.gaisr.pipeline/*slime-out*]
        (hunt-meta-gs "/data2/Bio/MetaG1/MoStos/RF00167/RF00167.sto.cm"))))
(future-done? *metag1-cmsearch-RF00167*)


(def *metag1-cmsearch-RF00174*
     (future
      (hunt-meta-gs "/data2/Bio/MetaG1/MoStos/RF00174/RF00174.sto.cm")))
(future-done? *metag1-cmsearch-RF00174*)


(def *metag1-cmsearch-RF00634*
     (future
      (hunt-meta-gs "/data2/Bio/MetaG1/MoStos/RF00634/RF00634.sto.cm")))
(future-done? *metag1-cmsearch-RF00634*)


(defmacro meta-hunt [name]
  `(def ~(symbol (str "*metag1-cmsearch-" name "*"))
        (future
         (hunt-meta-gs
          ~(str "/data2/Bio/MetaG1/MoStos/" name "/" name ".sto.cm")))))


(map #(list 'meta-hunt %)
     (drop 6 (sort (fs/listdir "/data2/Bio/MetaG1/MoStos"))))

(meta-hunt "RF01051")
(future-done? *metag1-cmsearch-RF01051*)

(meta-hunt "RF01054")
(future-done? *metag1-cmsearch-RF01054*)

(meta-hunt "RF01055")
(future-done? *metag1-cmsearch-RF01055*)

(meta-hunt "RF01692")
(future-done? *metag1-cmsearch-RF01692*)

(meta-hunt "RF01831")
(future-done? *metag1-cmsearch-RF01831*)




(def *base-fna-dir* "/data2/Bio/MetaG1/FastaFiles")

(def *metag1-fnas*
     (sort (fs/listdir *base-fna-dir*)))

(defn hunt-meta-gs [cms]
  (doseq [fna-pair (partition 2 2 nil *metag1-fnas*)]
    (let [start-time (. System (nanoTime))]
      (log> "metag1-runs" :info
            "MG1-Hunt: ~A , ~A"
            fna-pair
            ;; Make this part of logging so that it FORCES the total
            ;; completion of the computation of a pair.  This is
            ;; because log printer (IO) needs the result value to print
            ;; for each call to it.
            (cms&hitfnas->cmsearch-out
             cms (map #(fs/join *base-fna-dir* %) fna-pair)
             :par false :eval 1.0))
      (let [stop-time (. System (nanoTime))]
        (log> "metag1-runs" :info
              "         ~A secs, ~A mins"
              (/ (/ (double (- stop-time start-time)) 1000000.0)
                 1000.0)
              (/ (/ (double (- stop-time start-time)) 1000000.0)
                 60000.0))))))