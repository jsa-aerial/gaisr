
(in-ns 'edu.bc.bio.apps.qin)





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





;;; Various ROC data (sim, ecoli, bs) ---------- ;;;


;;(->> (fs/glob "/home/kaila/Bio/Test/ROC-data/Sim/RF*")(map fs/basename)sort)

(def our-sim-rnas
     ["RF00050" "RF00059" "RF00114" "RF00162" "RF00167"
      "RF00168" "RF00234" "RF00379" "RF00380" "RF00442"
      "RF00504" "RF00506" "RF00521" "RF00522" "RF00555"
      "RF00559" "RF00634" "RF01727" "RF01767" "RF01831"])

(def eval-sim-roc
     (eval-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Sim"
      :RNAs our-sim-rnas))

(def sccs-sim-roc
     (sccs-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Sim"
      :RNAs our-sim-rnas))


(def eval-ecoli-roc
     (eval-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Ecoli"
      :RNAs ["L1" "L10" "L20" "L4" "S1" "S15" "S2" "S4" "S7" "S8"]))

(def sccs-ecoli-roc
     (sccs-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Ecoli"
      :RNAs ["L1" "L10" "L20" "L4" "S1" "S15" "S2" "S4" "S7" "S8"]))


(def eval-bsub-roc
     (eval-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Bsub"
      :RNAs ["L20" "S15" "S4" "S6S"]))

(def sccs-bsub-roc
     (sccs-roc-data
      :base "/home/kaila/Bio/Test/ROC-data/Bsub"
      :RNAs ["L20" "S15" "S4" "S6S"]))



;;; Scores for ROC data ----------------------;;;

(def eval-sim-scores
     (->> eval-sim-roc
          (map #(ROC-opt-classifier-pt
                 % (fn[[nm cpt M]] (M :TPR)) (fn[[nm cpt M]] (M :FPR))))
          (sort-by first)))

(def sccs-sim-scores
     (->> sccs-sim-roc
          (map #(ROC-opt-classifier-pt
                 % (fn[[nm cpt M]] (M :TPR)) (fn[[nm cpt M]] (M :FPR))))
          (sort-by first)))


(def eval-ecoli-scores
     (->> eval-ecoli-roc
          (map #(ROC-opt-classifier-pt
                 % (fn[[nm cpt M]] (M :TPR)) (fn[[nm cpt M]] (M :FPR))))
          (sort-by first)))

(def sccs-ecoli-scores
     (->> sccs-ecoli-roc
          (map #(ROC-opt-classifier-pt
                 % (fn[[nm cpt M]] (M :TPR)) (fn[[nm cpt M]] (M :FPR))))
          (sort-by first)))


(def eval-bsub-scores
     (->> eval-bsub-roc
          (map #(ROC-opt-classifier-pt
                 % (fn[[nm cpt M]] (M :TPR)) (fn[[nm cpt M]] (M :FPR))))
          (sort-by first)))

(def sccs-bsub-scores
     (->> sccs-bsub-roc
          (map #(ROC-opt-classifier-pt
                 % (fn[[nm cpt M]] (M :TPR)) (fn[[nm cpt M]] (M :FPR))))
          (sort-by first)))



;;; Sorting out Sim data for final versions
;;;
(->> (map (fn[[sc [n & r]] [esc [en & er]]]
            (if (not= n en)
              (raise :type :name-mismatch :n n :en en)
              [n sc esc]))
          (sort-by #(-> % second first) sccs-sim-scores)
          (sort-by #(-> % second first) eval-sim-scores))
     (filter (fn[[n ssc esc]] (> esc ssc)))
     (sort-by second))


;;; Generate CSV ROC data for R ggplot
;;;
(->> (map (fn[[sc [n & r]] [esc [en & er]]]
            (if (not= n en)
              (raise :type :name-mismatch :n n :en en)
              [n sc esc]))
          (sort-by #(-> % second first) sccs-sim-scores)
          (sort-by #(-> % second first) eval-sim-scores))
     (filter (fn[[n ssc esc]] (>= esc ssc)))
     (sort-by third)
     (drop 3)
     (map first)
     #_(build-roc-curves :base "/home/kaila/Bio/Test/ROC-data/Sim" :RNAs)
     (csv-roc-data :base "/home/kaila/Bio/Test/ROC-data/Sim" :RNAs))



;;; Opt Cut Point generations

(defn roc->auc [roc-data]
  (let [rna (ffirst roc-data)
        data-maps (map third roc-data)
        pt-map (reduce (fn [M m]
                         (let [[tpr fpr] [(m :TPR) (m :FPR)]]
                           (cond (and (< fpr 1.0)
                                      (< (M fpr -1) tpr))
                                 (assoc M fpr tpr)

                                 (and (== fpr 1.0)
                                      (> (M fpr 2.0) tpr))
                                 (assoc M fpr tpr)

                                 :else M)))
                       {0.0 0.0, 1.0 1.0} data-maps)]
    [rna (->> pt-map (sort-by key)
              (#(interleave % (drop 1 %)))
              (partition-all 2)
              AUC ;;)]))
              (#(/ (- % 0.5) 0.5)))]))


(map roc->auc sccs-ecoli-roc)
(map roc->auc eval-ecoli-roc)

(map roc->auc sccs-bsub-roc)
(map roc->auc eval-bsub-roc)

(map roc->auc sccs-sim-roc)
(map roc->auc eval-sim-roc)



;;; Build tables
;;;

(defn rndto [flt dp]
  (let [d (expt 10.0 dp)]
    (-> flt (* d) round (/ d))))

(defn rnd24 [flt]
  (rndto flt 4))

(defn rnd25 [flt]
  (rndto flt 5))



(def multi-ctx
     {"S6S B.sub" true,
      "L1 E.coli" true,
      "RF00114" true,
      "RF00379" true,
      "RF00380" true,
      "RF00442" true,
      "RF00521" true,
      "RF00555" true,
      "RF00634" true})

(defn rna-nm [rna tx]
  (let [rna (if (not= tx "") (str rna " " tx) rna)]
    (if (multi-ctx rna) (str rna "\\$^M\\$") rna)))


;;; AUC table

(def auc-head
     "ncRNA & SCCS & Eval \\\\\\toprule")

(def auc-line
     "/rna/  & /sccs/ & /eval/ \\\\")

(def table-foot
     (str/join
      "\n" ["\\botrule" "\\end{tabular}}"]))

(defn gen-auc-line [tx sccs-data eval-data]
  (let [rows (map flatten
		  (transpose (map roc->auc sccs-data)
			     (map roc->auc eval-data)))]
    (doseq [[rna sccs _ eval] rows]
      (if (not= rna _)
	(raise :type :name-mismatch :args [rna sccs _ eval])
	(println (reduce (fn [L [re stg]]
			   (str/replace-re re stg L))
			 auc-line
			 [[#"/rna/" (rna-nm rna tx)]
			  [#"/sccs/" (-> sccs rnd24 str)]
			  [#"/eval/" (-> eval rnd24 str)]]))))))



;;; Main table

(def table-head
     "ncRNA  & Optctpt & TPR  & FPR  & PPV  & ACC & MCC \\\\\\toprule")

(def table-line
     "/rna/  & /ctpt/  & /tpr/  & /fpr/  & /ppv/  & /acc/  & /mcc/ \\\\")


(def bioinfo-table-template
     ["{\\scriptsize \\begin{tabular}{ll|lllll}\\toprule"
      "ncRNA  & Optctpt & TPR    & FPR    & PPV    & ACC    & MCC " "\\toprule"
      "rna    & /ctpt/  & /auc/  & /tpr/  & /spc/  & /ppv/  & /mcc/"
      "\\botrule"
      "\\end{tabular}}"])


(defn gen-line [tx data & {:keys [rndcpt?] :or {rndcpt? true}}]
  (let [[sc [rna cpt M]] data
        cpt (-> cpt Float.)
        cpt (if rndcpt? (-> cpt rnd24 str) (str cpt))
        TPR (-> M :TPR rnd24 str)
        FPR (-> M :FPR rnd24 str)
        PPV (-> M :ppv rnd24 str)
        ACC (-> M :ACC rnd24 str)
        MCC (-> M :MCC rnd24 str)
        rna (rna-nm rna tx)]
    (reduce (fn[L [re stg]]
              (str/replace-re re stg L))
            table-line
            [[#"/rna/" rna]
             [#"/ctpt/" cpt]
             [#"/tpr/" TPR]
             [#"/fpr/" FPR]
             [#"/ppv/" PPV]
             [#"/acc/" ACC]
             [#"/mcc/" MCC]])))

;;; SCCS

(doseq [[tx data] [["B.sub" (sort-by #(-> % second first) sccs-bsub-scores)]
                   ["E.coli" (sort-by #(-> % second first) sccs-ecoli-scores)]
                   ["" (sort-by #(-> % second first) sccs-sim-scores)]]]
  (doseq [d data]
    (println (gen-line tx d))))

;;; Eval

(doseq [[tx data] [["B.sub" (sort-by #(-> % second first) eval-bsub-scores)]
                   ["E.coli" (sort-by #(-> % second first) eval-ecoli-scores)]
                   ["" (sort-by #(-> % second first) eval-sim-scores)]]]
  (doseq [d data]
    (println (gen-line tx d :rndcpt? false))))



;;; Other pictures / charts

