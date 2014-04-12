
(in-ns 'edu.bc.bio.apps.qin)




(->> (map (fn[[sc [n & r]] [esc [en & er]]]
            (if (not= n en)
              (raise :type :name-mismatch :n n :en en)
              [n sc esc]))
          (sort-by #(-> % second first) sccs-scores)
          (sort-by #(-> % second first) eval-scores))
     (filter (fn[[n ssc esc]] (> esc ssc)))
     (sort-by second))

(->> (map (fn[[sc [n & r]] [esc [en & er]]]
            (if (not= n en)
              (raise :type :name-mismatch :n n :en en)
              [n sc esc]))
          (sort-by #(-> % second first) sccs-scores)
          (sort-by #(-> % second first) eval-scores))
     (filter (fn[[n ssc esc]] (>= esc ssc)))
     (sort-by third)
     (drop 3)
     (map first)
     #_(build-roc-curves :base "/home/kaila/Bio/Test/ROC-data/Sim" :RNAs)
     (csv-roc-data :base "/home/kaila/Bio/Test/ROC-data/Sim" :RNAs))





(let [base "/home/kaila/Bio/Test/ROC-data/Ecoli/S8"
      final (->> (fs/join base "CSV/*-final.ent") fs/glob first)
      good (->> final slurp csv/parse-csv
                (take-while (fn[[n jsd]] (< (Double. jsd) 0.919))))
      bad (->> final slurp csv/parse-csv
               (drop-while (fn[[n jsd]] (< (Double. jsd) 0.919))))
      good-file (fs/join base "TotPosNeg/good-finals.ent")
      bad-file (fs/join base "TotPosNeg/bad-finals.ent")]
  (io/with-out-writer good-file
    (doseq [[n jsd] good]
      (println n))
    (println))
  (io/with-out-writer bad-file
    (doseq [[n jsd] bad]
      (println n))
    (println)))




;;; One Sim
(time (doall (build-sccs-pos-neg-sets
              (->> (fs/directory-files
                    "/home/kaila/Bio/Test/ROC-data/Sim/RF00050/CSV"
                    ".cmsearch.csv") first))))
(build-eval-pos-neg-sets
 "/home/kaila/Bio/Test/ROC-data/Sim/RF00050/CSV/RF00050-seed-NC.sto.all-sim-db.cmsearch.csv")


(binding [sccs-cut-points (->> (range -0.25 4.0 0.1) (map logistic) (take 25))]
  (doseq [rf ["RF00379"]] ;; "RF00380" "RF00442" "RF00634"]]
    (time (doall (build-sccs-pos-neg-sets
                  (->> (fs/directory-files
                        (fs/join "/home/kaila/Bio/Test/ROC-data/Sim" rf "CSV" )
                        ".cmsearch.csv")
                       first))))))

(doseq [rf ["RF00522" "RF01767" "RF01831"]]
  (time (doall (build-sccs-pos-neg-sets
                (->> (fs/directory-files
                      (fs/join "/home/kaila/Bio/Test/ROC-data/Sim" rf "CSV" )
                      ".cmsearch.csv")
                     first)))))



(defn csv-bkup->ent [csv-bkup & {:keys [entfile]}]
  (let [csv (->> csv-bkup slurp
                 ((fn[x] (spit "/tmp/x.csv" x) "/tmp/x.csv")))
        entfile (if entfile
                  (fs/fullpath entfile)
                  (fs/replace-type csv-bkup ".ent"))]
    (io/with-out-writer entfile
      (doseq [ent (get-entries csv)]
        (println ent)))))

(let [base "/home/kaila/Bio/Test/ROC-data/Sim"]
  (doseq [rf (sort (->> (fs/join base "RF*") fs/glob (map fs/basename)))]
    (let [TotDir (fs/join base rf "TotPosNeg")]
      (when (not (fs/exists? TotDir)) (fs/mkdir TotDir))
      (csv-bkup->ent
       (->> (fs/join base rf "/CSV/*bkup*") fs/glob first slurp
            ((fn[x] (spit "/tmp/x.csv" x) "/tmp/x.csv")))
       :entfile (fs/join TotDir (str rf "-neg.ent")))
      (sto->ent
       (fs/join base rf (str rf "-seed-NC.sto"))
       :entfile (fs/join TotDir (str rf "-pos.ent"))))))


;;; Old, New
(let [all (->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")
               (map #(->> % fs/basename (str/split #"-") first)))
      new (->> (fs/glob "/data2/BioData/RFAM-NC/NC/*.sto")
               (map #(->> % fs/basename (str/split #"-") first)))
      old (set/difference (set all) (set new))]
  (sort old))




(for [rf ["RF00380" "RF00442" "RF00634"]
      :let [base (fs/join "/home/kaila/Bio/Test/ROC-data/Sim" rf)
            stos (fs/glob (fs/join base "CLU/*.sto"))
            chart-dir (fs/join base "CLU/Charts")]]
  (doseq [sto stos]
    (sccs/save-hit-context-delta sto :plot chart-dir)))







(->> (sccs-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Ecoli"
      :RNAs ["L1" "L10" "L20" "L4" "S1" "S15" "S2" "S4" "S7" "S8"])
     (map #(map third %))
     (map #(reduce (fn [M m]
		     (let [[tpr fpr] [(m :TPR) (m :FPR)]]
		       (if (< (M fpr -1) tpr)
			 (assoc M fpr tpr)
			 M)))
		   {0.0 0.0, 1.0 1.0} %))
     (map #(sort-by key %))
     (map #(->> (interleave % (drop 1 %)) (partition-all 2)))
     (map AUC))