;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                    S E Q U T I L S . D I S T S                           ;;
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

(ns edu.bc.bio.sequtils.dists

  "Various sequence, sequence families, alignments and alignment
   families distributions.  Intended for use as \"background\"
   distributions in relative entropy and other such comparaitive or
   measure calculations.  These distributions are focused on microbial
   NC genomes and RFAM alignment families.  At least to start.
   Others, such as metagenomic results and NZ genomes could be here in
   future."

  (:require [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.io :as io]
            [incanter.core]
            [incanter.charts]
            [edu.bc.fs :as fs]
            [edu.bc.utils.graphs :as gr]
            [edu.bc.utils.clustering :as clu])
  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.info-theory
        [clojure.pprint
         :only [cl-format]]
        ))




;;; ----------------------------------------------------------------------
;;;
;;; Some bacground frequency and probability distributions for all NC
;;; genomes and RFAM families.  May well be better to place in
;;; database and be pulled in and initialized at startup time...
;;;

(defn print-dist
  "Prints a distribution DIST, a frequency or probability distribution
   map, in sorted order of most occurances / highest probability to
   lowest.  One item per line.
  "
  [dist]
  (doseq [[k v] (sort-by second > dist)]
    (print " ") (pr k) (print " ") (pr v) (println \,)))


(defn prn-dist
  "Print readably the distribution DIST to file FILE.  This means that
   FILE can be used as a secondary cache and so can be simply read
   with Clojures reader to get the internal representation back."
  [dist file]
  (io/with-out-writer (fs/fullpath file)
    (prn dist)))



;;; (bg-freqs-probs 2 "/data2/BioData/GenomeSeqs"
;;;                 :ftypes [".fna"] :norm false :par 4)

(def +bg-nc-genomes+
     ^{:doc
       "Background dinucleotide Frequencies for all NC microbial
        genomes.  In addition to the 16 \"standard\" 4 base
        pair (A,T,G,C) permutations, includes many non standard
        pairs (see +NONSTD-RNA+"}

     {"GC" 395069986,
      "CG" 377952477,
      "TT" 326305148,
      "AA" 326263282,
      "CC" 306796217,
      "GG" 306596839,
      "AT" 299033486,
      "CA" 291400632,
      "TG" 291273508,
      "TC" 287515297,
      "GA" 287493041,
      "AG" 251060051,
      "CT" 251035265,
      "AC" 237803474,
      "GT" 237723100,
      "TA" 209003358,
      "NN" 11849, "YT" 693, "NA" 688, "AN" 676, "TN" 664, "NT" 626, "NG" 480,
      "NC" 470, "GN" 468, "CN" 458, "AR" 403, "TY" 359, "CY" 335, "GY" 313,
      "RA" 295, "RG" 245, "RC" 229, "YG" 227, "AY" 227, "CR" 222, "CS" 189,
      "MC" 184, "CM" 177, "YC" 177, "SC" 172, "AK" 171, "RT" 170, "WA" 162,
      "TW" 156, "AM" 155, "GS" 155, "GR" 153, "TR" 150, "KG" 149, "YA" 143,
      "SG" 143, "MA" 140, "AS" 139, "GK" 135, "KA" 133, "TK" 127, "KC" 126,
      "ST" 125, "KT" 121, "MT" 110, "CK" 108, "WT" 107, "SA" 106, "GM" 102,
      "AW" 101, "CW" 99, "MG" 99, "TM" 93, "WG" 84, "TS" 76, "WC" 74, "GW" 68,
      "RR" 25, "MM" 18, "SS" 17, "KK" 14, "YY" 14, "YW" 14, "KR" 13, "SM" 12,
      "WW" 11, "SY" 11, "MY" 11, "WR" 11, "RK" 10, "WS" 10, "WY" 9, "RY" 9,
      "KY" 9, "YK" 9, "YR" 9, "SR" 9, "KM" 8, "SW" 8, "KS" 8, "MW" 8, "RM" 8,
      "BG" 7, "CH" 7, "MR" 7, "RW" 7, "GB" 7, "RS" 6, "VA" 6, "BT" 6, "CV" 6,
      "WM" 6, "SK" 6, "YS" 6, "CB" 6, "AD" 5, "MS" 5, "KW" 5, "VC" 5, "AV" 5,
      "HC" 5, "CD" 4, "TV" 4, "DG" 4, "HT" 4, "TB" 4, "GV" 4, "VG" 4, "YM" 4,
      "BA" 4, "DC" 4, "HH" 3, "GH" 3, "WK" 3, "SN" 3, "VT" 3, "AB" 2, "DT" 2,
      "TD" 2, "YN" 2, "DA" 2, "GD" 2, "NK" 2, "MK" 2, "HG" 2, "NM" 2, "BC" 1,
      "MN" 1, "VW" 1, "KN" 1, "VY" 1, "NR" 1, "AH" 1, "DM" 1, "MV" 1, "NW" 1,
      "HR" 1, "NY" 1, "YH" 1, "BR" 1, "TH" 1, "WN" 1, "HA" 1, "RN" 1, "VS" 1,
      "WV" 1})

(def +bg-nc-genomes-probs+
     ^{:doc
       "Background dinucleotide probabilities for all NC microbial
       genomes."}
     (probs +bg-nc-genomes+))


(def +bg-nc-genomes-sym+
     ^{:doc
       "Background dinucleotide Frequencies for all NC microbial
        genomes where pairs are treated symmetrically, i.e., XY is
        considerted = YX.  In addition to the 10 \"standard\" 4 base
        pair (A,T,G,C) symmetric permutations, includes many non
        standard pairs (see +NONSTD-RNA+"}

     {"GC" 773022463,
      "GA" 538553092,
      "TC" 538550562,
      "CA" 529204106,
      "TG" 528996608,
      "AT" 508036844,
      "TT" 326305148,
      "AA" 326263282,
      "CC" 306796217,
      "GG" 306596839,
      "NN" 11849, "NA" 1364, "TN" 1290, "YT" 1052, "NG" 948, "NC" 928,
      "AR" 698, "GY" 540, "CY" 512, "RC" 451, "RG" 398, "AY" 370, "CS" 361,
      "MC" 361, "RT" 320, "AK" 304, "GS" 298, "AM" 295, "KG" 284, "TW" 263,
      "WA" 263, "TK" 248, "AS" 245, "KC" 234, "MT" 203, "ST" 201, "GM" 201,
      "CW" 173, "WG" 152, "RR" 25, "KR" 23, "YW" 23, "MM" 18, "RY" 18,
      "KY" 18, "WR" 18, "WS" 18, "SS" 17, "SY" 17, "SM" 17, "MY" 15, "RM" 15,
      "SR" 15, "KK" 14, "YY" 14, "BG" 14, "KS" 14, "MW" 14, "CH" 12, "WW" 11,
      "VA" 11, "CV" 11, "KM" 10, "BT" 10, "CD" 8, "KW" 8, "GV" 8, "TV" 7,
      "AD" 7, "CB" 7, "DG" 6, "BA" 6, "GH" 5, "HT" 5, "DT" 4, "HH" 3,
      "YN" 3, "SN" 3, "NK" 3, "NM" 3, "VW" 2, "NR" 2, "AH" 2, "NW" 2,
      "VY" 1, "DM" 1, "MV" 1, "HR" 1, "YH" 1, "BR" 1,
      "VS" 1})

(def +bg-nc-genomes-sym-probs+
     ^{:doc
       "Background dinucleotide probabilities for symmetric
       frequencies for all NC microbial genomes."}
     (probs +bg-nc-genomes-sym+))


;;; (bg-freqs-probs 1 "/data2/Bio/RFAM" :ftypes [".sto"] :norm true :par 2)
(def +bg-freqs-probs-all-rfam+
     (let [fs { "." 634856, "A" 157097, "G" 154623, "C" 130024, "U" 128626,
                "N" 7, "Y" 1, "S" 1}]
       [fs (probs fs)]))

;;; (bg-freqs-probs 2 "/data2/Bio/RFAM" :ftypes [".sto"] :norm true :par 2)
(def +bg-freqs-probs-all-rfam-di+
     (let [fs { ".." 570773, "AA" 46303, "GG" 40160, "GA" 36784, "AG" 34170,
                "GC" 32701, "UU" 30849, "CC" 30216, "UG" 29804, "CG" 29361,
                "GU" 28772, "CA" 28589, "AU" 28537, "AC" 27831, "UC" 27120,
                "UA" 25212, "CU" 25150, ".G" 20125, "A." 18870, ".A" 18021,
                "C." 15578, "G." 15092, "U." 14374, ".U" 14005, ".C" 11403,
                "CN" 3, "N." 2, "NU" 2, "NA" 2, "UN" 2, "SU" 1, "UY" 1,
                "GN" 1, "AN" 1, "NC" 1, "YU" 1, "US" 1}]
       [fs (probs fs)]))

;;; (bg-freqs-probs 2 "/data2/Bio/RFAM"
;;;                 :fsps-fn cc-combins-freqs-probs
;;;                 :cols true :ftypes [".sto"] :norm true :par 2)
(def +bg-freqs-probs-all-rfam-cols+
     (let [fs {".." 165988848, "GG" 15598223, "AA" 14862080, "CC" 11859618,
               "UU" 8978659, "AG" 6619716, "UA" 6160495, "UC" 5644286,
               "A." 4008987, ".U" 3875698, "AC" 3812954, "CG" 3769734,
               "GU" 3742386, "G." 3532558, "C." 3352204, "GA" 340019,
               "CU" 315747, ".A" 288295, ".G" 243124, ".C" 142407,
               "GC" 72061, "U." 59768, "CA" 51845, "AU" 30890,
               "UG" 14227, "N." 353, "AN" 244, "CY" 178, "GN" 117,
               "CN" 99, ".S" 85, "NU" 75, "GS" 10, "SU" 9, "SA" 6,
               "CS" 4, ".N" 3, "NC" 2, "UN" 2, "NG" 1}]
       [fs (probs fs)]))

;;; (bg-freqs-probs 2 "/data2/Bio/RFAM"
;;;                 :fsps-fn cc-combins-freqs-probs
;;;                 :cols true :ftypes [".sto"] :sym? true :norm true :par 2)
(def +bg-freqs-probs-all-rfam-cols-sym+
     (let [fs {".." 165988848, "GG" 15598223, "AA" 14862080, "CC" 11859618,
               "UU" 8978659, "AG" 6959735, "UA" 6191385, "UC" 5960033,
               "A." 4297282, ".U" 3935466, "AC" 3864799, "CG" 3841795,
               "G." 3775682, "GU" 3756613, "C." 3494611, "N." 356,
               "AN" 244, "CY" 178, "GN" 118, "CN" 101, ".S" 85, "NU" 77,
               "GS" 10, "SU" 9, "SA" 6, "CS" 4}]
       [fs (probs fs)]))



(def
 ^{:doc
   "Memoized frequency and probability distributions over RFAM
    families.  The results returned are maps from an RFAM family
    identified by its string name (e.g., \"RF00050\") to pairs of the
    families frequency and probability distributions (as maps).  Keys
    are determined by n: 1 = mononucleotide, 2 = dinucleotide, etc.
    All nucleotide keys are strings.

    Ex:
    (sort (keys (bg-freqs-probs-rfam-by-family 1)))
    => ; A second or so
    (\"RF00050\" \"RF00057\" \"RF00059\" \"RF00114\" \"RF00140\" \"RF00162\"
     \"RF00167\" \"RF00168\" \"RF00174\" \"RF00234\" \"RF00379\" \"RF00380\"
     \"RF00442\" \"RF00504\" \"RF00506\" \"RF00514\" \"RF00515\" \"RF00521\"
     \"RF00522\" \"RF00555\" \"RF00556\" \"RF00557\" \"RF00558\" \"RF00559\"
     \"RF00634\" \"RF01051\" \"RF01054\" \"RF01055\" \"RF01057\" \"RF01070\"
     \"RF01385\" \"RF01402\" \"RF01482\" \"RF01510\" \"RF01689\" \"RF01692\"
     \"RF01693\" \"RF01694\" \"RF01704\" \"RF01725\" \"RF01727\" \"RF01739\"
     \"RF01745\" \"RF01767\" \"RF01769\" \"RF01793\" \"RF01826\" \"RF01831\")

     ((bg-freqs-probs-rfam-by-family 1) \"RF00050\")
     => ; Second time - cached and immediate
     [{\"C\" 4398, \"A\" 5034, \"G\" 6278, \"U\" 4212, \".\" 12344}
      {\".\" 0.3825698878075993, \"U\" 0.1305398871877518,
       \"G\" 0.1945701357466063, \"A\" 0.1560156201574413,
       \"C\" 0.1363044691006013}]

     (jensen-shannon
       (second ((bg-freqs-probs-rfam-by-family 2 :cols true) \"RF00050\"))
       (second +bg-freqs-probs-all-rfam-cols+))
     => 0.11288
   "

       :arglists '[n & {:keys [fsps-fn cols sym?]
                         :or {fsps-fn cc-freqs-probs cols false sym? false}}]}

     bg-freqs-probs-rfam-by-family

     (memoize
      (fn [n & {:keys [fsps-fn cols sym?]
                :or {fsps-fn cc-freqs-probs cols false sym? false}}]
        (if (and (= n 2)
                 (= fsps-fn cc-combins-freqs-probs)
                 cols)
          (io/with-in-reader "/data2/BioData/Archives/rfam-freqs-probs.clj"
            (read))
          (into
           {} (map (fn[file]
                     (let [nm (re-find #"^RF[0-9]+" (fs/basename file))
                           sqs (read-aln-seqs file :cols cols)
                           [fs ps] (take 2 (rest (seqs-freqs-probs
                                                  n sqs :fsps-fn fsps-fn
                                                  :nogaps false
                                                  :norm true :par 2)))]
                       [nm [fs ps]]))
                   (sort (fs/directory-files "/data2/Bio/RFAM" ".sto"))))))))


;;;(bg-freqs-probs-rfam-by-family 2)
;;;((bg-freqs-probs-rfam-by-family 2) "RF00050")
;;;
;;;
;;; This next one is quite expensive - several minutes and so we
;;; probably should secondarily cache this - that's where the DB
;;; should come in, but maybe just hack it to start with by writing it
;;; to a file...
;;;
;;;(bg-freqs-probs-rfam-by-family 2 :cols true :fsps-fn cc-combins-freqs-probs)

;;;(prn-dist
;;; (bg-freqs-probs-rfam-by-family 2 :cols true :fsps-fn cc-combins-freqs-probs)
;;; "/data2/BioData/Archives/rfam-freqs-probs.clj")

;;; (map #(do [% (jensen-shannon (second ((bg-freqs-probs-rfam-by-family
;;;                                        2 :cols true
;;;                                        :fsps-fn cc-combins-freqs-probs) %))
;;;                              (second +bg-freqs-probs-all-rfam-cols+))])
;;;      (sort (keys (bg-freqs-probs-rfam-by-family 1))))

;;; (def +ecr-rfam-js+
;;;      (for [sto (fs/directory-files "/data2/Bio/ECRibLeaders/MoStos" ".sto")
;;;            :let [nm (first (str/split #"\." (fs/basename sto)))
;;;                  stodist (second (bg-freqs-probs
;;;                                   2 [sto] :cols true
;;;                                   :fsps-fn cc-combins-freqs-probs
;;;                                   :par 2))]]
;;;        [nm (sort-by
;;;             second
;;;             (map #(do [% (jensen-shannon
;;;                           stodist
;;;                           (second ((bg-freqs-probs-rfam-by-family
;;;                                     2 :cols true
;;;                                     :fsps-fn cc-combins-freqs-probs) %)))])
;;;                  (sort (keys (bg-freqs-probs-rfam-by-family 1)))))]))

;;; (time
;;;  (bg-freqs-probs
;;;   2 (take 1 (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa"))
;;;   :par 4))


#_(time
   (io/with-out-writer "/data2/Bio/MetaG1/lHmax.clj"
     (doseq [mg-fa (take 20 (sort (fs/directory-files
                                   "/data2/Bio/MetaG1/FastaFiles" ".fa")))]
       (prn [(fs/basename mg-fa)
             (map #(let [[nm sq] %]
                     [nm (count sq) (-> sq count log4)
                      (for [l (range 1 20)]
                        [l (count (keys (reduce
                                         (fn[m [k v]]
                                           (if (> v 1)
                                             (assoc m k v) m))
                                         {} (freqn l sq))))])])
                  (read-seqs mg-fa :info :both))]))))



;;;--------- REFACTOR - The following needs to go into info-theory or ??? ;;;

(defn expected-qdict [q-1 q-2 & {:keys [alpha] :or {alpha (alphabet :rna)}}]
  (reduce (fn[m lmer]
            (let [l (dec (count lmer))
                  x (subs lmer 1)
                  y (subs lmer 0 l)
                  z (subs lmer 1 l)]
              (if (and (q-1 x) (q-1 y) (q-2 z))
                (assoc m lmer (/ (* (q-1 x) (q-1 y)) (q-2 z)))
                m)))
          {} (for [x (keys q-1) a alpha] (str x a))))


(defn freq-xdict-dict
  [q sq]
  (let [ext-sq (str sq "X")
        residue-key (subs ext-sq (- (count sq) (- q 1)))
        q-xfreq-dict (freqn q ext-sq)
        q-freq-dict (dissoc q-xfreq-dict residue-key)]
    [q-xfreq-dict q-freq-dict]))

(defn q-1-dict
  ([q-xdict]
     (reduce (fn[m [k v]]
               (let [l (count k)
                     pre (subs k 0 (dec l))]
                 (assoc m pre (+ (get m pre 0) v))))
             {} q-xdict))
  ([q sq]
     (probs (dec q) sq)))

(defn q1-xdict-dict
  [q sq & {:keys [ffn] :or {ffn probs}}]
  (let [[q-xfreq-dict q-freq-dict] (freq-xdict-dict q sq)
        q-xpdf-dict (probs q-xfreq-dict)
        q-pdf-dict (probs q-freq-dict)]
    {:xfreq q-xfreq-dict :xpdf q-xpdf-dict
     :freq q-freq-dict :pdf q-pdf-dict}))



(defn reconstruct-dict
  [l sq & {:keys [alpha] :or {alpha (alphabet :rna)}}]
  {:pre [(> l 2)]}
  (let [q (dec l)
        qmaps (q1-xdict-dict q sq)
        [q-xdict q-dict] (map qmaps [:xpdf :pdf])
        q-1dict (q-1-dict q-xdict)]
    (expected-qdict q-dict q-1dict :alpha alpha)))


(defn max-qdict-entropy
  [q & {:keys [alpha] :or {alpha (alphabet :rna)}}]
  (let [base (count alpha)]
    (* q (log2 base))))

(defn informativity
  ([q sq]
     (- (max-qdict-entropy q) (entropy (probs q sq))))
  ([q-dict]
     (let [q (count (first (keys q-dict)))]
       (- (max-qdict-entropy q) (entropy q-dict)))))


(defn limit-entropy
  [q|q-dict sq|q-1dict &
   {:keys [alpha NA] :or {alpha (alphabet :rna) NA -1.0}}]
  {:pre [(or (and (integer? q|q-dict)
                  (or (string? sq|q-1dict) (coll? sq|q-1dict)))
             (and (map? q|q-dict) (map? sq|q-1dict)))]}

  (if (map? q|q-dict)
    (let [q-dict q|q-dict
          q-1dict sq|q-1dict]
      (/ (- (entropy q-dict) (entropy q-1dict))
         (log2 (count alpha))))

    (let [q q|q-dict
          sq sq|q-1dict
          lgcnt (log2 (count alpha))]
      (if (= q 1)
        (/ (entropy (probs 1 sq)) lgcnt)

        (let [qmaps (q1-xdict-dict q sq)
              [q-xdict q-dict] (map qmaps [:xpdf :pdf])
              q-1dict (q-1-dict q-xdict)]
          (if (< (count q-dict) (count q-1dict))
            (if (fn? NA) (NA q-dict q-1dict) NA)
            (/ (- (entropy q-dict) (entropy q-1dict)) lgcnt)))))))


(defn limit-informativity
  ([q sq]
     )
  ([q-dict]
     ))


(defn CREl
  [l sq & {:keys [limit alpha]
           :or {limit 15 alpha (alphabet :rna)}}]
  {:pre [(> l 2)]}
  (sum (fn[k]
         (catch-all (DX||Y
                     (probs k sq)
                     (reconstruct-dict k sq :alpha alpha))))
       (range l (inc limit))))


(defn information-capacity
  [q sq & {:keys [cmpfn] :or {cmpfn jensen-shannon}}]
  (catch-all (cmpfn (probs q sq)
                    (reconstruct-dict q sq))))


(defn hybrid-dictionary
  [l sqs]
  {:pre [(or (string? sqs) (coll? sqs))]}
  (let [sqs (if (-> sqs first map?)
              sqs
              (degap-seqs (if (coll? sqs) sqs (read-seqs sqs))))
        cnt (count sqs)
        par (max (floor (/ cnt 10)) 2)
        dicts (if (-> sqs first map?) sqs (xfold #(probs l %) sqs))
        hybrid (apply merge-with +
                      (xfold (fn[subset] (apply merge-with + subset))
                             (partition-all (/ (count dicts) par) dicts)))]
    (reduce (fn[m [k v]] (assoc m k (double (/ v cnt))))
            {} hybrid)))


;;; (pxmap (fn[subset] (apply merge-with + subset))
;;;                              par (partition-all
;;;                                   (/ (count dicts) par)
;;;                                   dicts))
;;; (xfold (fn[subset] (apply merge-with + subset))
;;;                          (partition-all (/ (count dicts) par) dicts))

;;; 1774444 the number of keys!!

;;; (/ 186112 2)
;;; 57 253 960 bases
;;;
;;; 4 683 023 485 bases
;;; (/ 7412 2) => 3706
;;;
;;; (/ 4683023485 57253960.0)
;;; (/ 7512557802 57253960.0)


(defn ctx-seq
  [entry & {:keys [directed ldelta rdelta delta ddel] :or {directed true}}]
  {:pre [(or delta (and (not directed) (or ldelta rdelta)))]}
  (let [ldelta (or delta ldelta 0)
        rdelta (or delta rdelta 0)]
    (if (not directed)
      (gen-name-seq entry :ldelta ldelta :rdelta rdelta)
      (let [ddel (or ddel
                     (->> entry entry-parts second
                          (#(apply - %)) abs (#(/ % 4)) ceil))
            +? (= 1 (->> entry (pos \-) count))
            ldelta (if +? ddel delta)
            rdelta (if +? delta ddel)]
        (gen-name-seq entry :ldelta ldelta :rdelta rdelta)))))

(defn get-adjusted-seqs
  ""
  [entries delta & {:keys [directed ldelta rdelta ddel] :or {directed true}}]
  (let [ddel (or ddel (if (= delta 0) 0 nil))]
    (xfold #(ctx-seq % :directed directed
                     :delta delta :ddel ddel :ldelta ldelta :rdelta rdelta)
           entries)))

(defn cre-samples
  [seqs & {:keys [alpha xlate directed ldelta rdelta delta cnt limit]
           :or {alpha (alphabet :rna)
                directed true cnt 5 limit 14}}]
  {:pre [(or delta (and (not directed) ldelta rdelta))]}
  (let [entries (if (coll? seqs) seqs (get-entries seqs))
        cnt (min cnt (count entries))
        ;;ldelta (or delta ldelta)
        ;;rdelta (or delta rdelta)
        name-seq-pairs (get-adjusted-seqs entries delta
                                          :ldelta ldelta :rdelta rdelta)]
    (for [[nm sq] (random-subset name-seq-pairs cnt)
          :let [sq (if xlate (seqXlate sq :xmap xlate) sq)]]
      (xfold (fn[l] [l (CREl l sq :alpha alpha)
                     (count sq) nm])
             (range 3 (inc limit))))))


(defn render-chart
  ([location xs ys xlabel ylabel title & {:keys [series-label legend]}]
     (let [chart (if (or series-label legend)
                   (incanter.charts/scatter-plot
                    xs ys
                    :x-label xlabel :y-label ylabel :title title
                    :series-label series-label :legend legend)
                   (incanter.charts/scatter-plot
                    xs ys
                    :x-label xlabel :y-label ylabel :title title))]
       (render-chart chart location)))
  ([chart location]
     (if (= location :display)
       (incanter.core/view chart)
       (incanter.core/save chart location))))

(defn plot-cres
  [cre-samples & {:keys [file]}]
  (doseq [cres cre-samples]
    (let [ks (map first cres)
          cnt (third (first cres))
          cres (map second cres)]
      (render-chart
       (if file file :display)
       ks cres
       "feature length"
       "Commulative Relative Entropy"
       (str "CRE(l, F/F^), len: " cnt
            " Fmax: " (round (ln cnt)))))))


(defn res-word-size
  ""
  [seqs & {:keys [alpha xlate directed ldelta rdelta delta cnt limit crecut]
           :or {alpha (alphabet :rna)
                directed true cnt 5 limit 14 crecut 0.10}}]
  (let [cres (cre-samples seqs :delta delta
                          :xlate xlate :alpha alpha
                          :limit limit :cnt cnt)
        wz (reduce (fn [res v]
                     (/ (+ res
                           (some #(when (< (second %) crecut) (first %)) v))
                        2.0))
                   0.0 cres)]
    [(ceil wz) cres]))


(defn sto-re-dists
  ""
  [wz candidate-file sto-file
   & {:keys [refn delta order xlate par]
      :or {refn DX||Y delta 5000 order :up par 10}}]
  (let [comp (if (= order :up) < >)
        cfile (fs/fullpath candidate-file)
        ctype (fs/ftype cfile)
        prot-name (->> cfile
                       (str/split #"/") last (str/split #"\.") first
                       (str/replace-re #"_" " "))

        entries (get-entries cfile)
        sqs (->> entries
                 (#(get-adjusted-seqs % delta :ddel 0))
                 (map second)
                 (#(if xlate (seqXlate % :xmap xlate) %)))

        refsqs (->> sto-file
                    get-entries
                    (#(get-adjusted-seqs % delta :ddel 0))
                    (map second)
                    (#(if xlate (seqXlate % :xmap xlate) %)))

        sto-hd (hybrid-dictionary wz refsqs)
        dicts (xfold #(probs wz %) sqs)
        stods (xfold #(refn % sto-hd) dicts)]

    [(sort-by second comp (set (map (fn[en d] [en d]) entries stods)))
     prot-name
     (count sqs)]))


(defn nm-res-dist [nm-res]
  (let [frq (reduce (fn[m [nm re]]
                      (let [re (-> re (* 1000) (+ 0.1) round (/ 1000.0))]
                        (assoc m re (conj (get m re []) nm))))
                    {} nm-res)
        cnt (double (count nm-res))]
    (map (fn[[re nms]] [re (double (/ (count nms) cnt)) nms]) frq)))

(defn nm-res-cdf [nm-res]
  (let [re-dist (nm-res-dist nm-res)
        re-sorted (sort-by first re-dist)
        pre-sq (map (fn[[re p nms]] [p re]) re-sorted)]
    (fn[x]
      (reduce (fn[y [p re]]
                (if (<= re x) (+ y p) y))
              0.0 pre-sq))))

(defn re-cdf-cut [nm-res & {:keys [Dy Mre] :or {Dy 0.0 Mre 0.9}}]
  (let [res (sort (map first (nm-res-dist nm-res))) ; only the SET of res
        Fx (nm-res-cdf nm-res)
        pts (sort (map second nm-res))
        cdf-cut (+ Dy 0.5) ; 0.5 = median of Fx by definition
        re-cutoff (reduce (fn[v re]
                            (if (>= Dy 0.0)
                              ;; Straight CDF pick
                              (if (<= (Fx re) cdf-cut) re v)
                              ;; Add that re <= Max RE, to ensure late
                              ;; runs don't derail with CDF picking
                              ;; many bad, since nearly all are bad
                              (if (and (<= re Mre) (<= (Fx re) cdf-cut)) re v)))
                          0.0 res)
        ;;_ (print :re-cutoff1 re-cutoff :Fres (first res) "/ ")
        ;; A _single_ good score can occur an 'anomalous' number of times
        re-cutoff (let [fres (first res)]
                    (cond (not (zero? re-cutoff)) re-cutoff
                          (<= fres Mre) fres
                          :else 0.0))
        ;;_ (print :re-cutoff2 re-cutoff "/ ")
        ;; Round to nearest thousandth
        re-cutoff (/ (round (* 1000 (+ re-cutoff 0.001))) 1000.0)]
    #_(println :re-cutoff re-cutoff)
    (reduce (fn[cutpt re] (if (<= re re-cutoff) (inc cutpt) cutpt))
            0 (map second nm-res))))


(defn plot-re-pmf
  [nm-res  & {:keys [file]}]
  (let [pmf-info (map (fn[[re p v]] [re (count v)]) (nm-res-dist nm-res))
        pmf-dist (into {} pmf-info)
        xs (sort (map first pmf-info))]
    (render-chart
     (if file file :display)
     xs (map #(get pmf-dist % 0.0) xs)
     "RE value" "Sq/RE count" "Sq RE distribution"
     :series-label "sq/hbrid JSD PMF"
     :legend true)))

(defn plot-re-cdf
  [nm-res & {:keys [file]}]
  (let [Fx (nm-res-cdf nm-res)
        step 0.005
        xs (for [i (range 200)] (* i step))]
    (render-chart
     (if file file :display)
     xs (map Fx xs)
     "JSD RE" "Fre-dist(x)" "FFP CDF plot"
     :series-label "sq/hbrid JSD CDF"
     :legend true)))


(defn select-cutpoint
  [re-points & {:keys [area Dy Mre] :or {Mre 0.899}}]
  (if Dy
    (re-cdf-cut re-points :Dy Dy :Mre Mre)
    (let [s (* (sum re-points) area)]
      (first
       (reduce (fn[[x v] re]
                 (if (< v s) [(inc x) (+ v re)] [x v]))
               [0 0]
               re-points)))))

(defn get-good-candidates
  [nm-re-coll cutpoint & {:keys [ends] :or {ends :low}}]
  (let [total (count nm-re-coll)
        [good bad] (cond
                    (= ends :low)
                    [(take cutpoint nm-re-coll)
                     (drop cutpoint nm-re-coll)]

                    (= ends :high)
                    [(drop (- total cutpoint) nm-re-coll)
                     (take (- total cutpoint) nm-re-coll)]

                    :else ; :both
                    [(concat
                      (take cutpoint nm-re-coll)
                      (drop (- total cutpoint) nm-re-coll))
                     (take (- total cutpoint cutpoint)
                           (drop cutpoint nm-re-coll))])]
    [good bad]))


(defn selection-perf
  [nm-re-coll ent-file cutpoint & {:keys [ends good] :or {ends :low}}]
   (let [candidates (->> (io/read-lines ent-file)
                         (map entry-parts)
                         (map (fn[[nm [s e] sd]]
                                [(str nm "/" s "-" e "/" sd) 1])))
         total (count nm-re-coll)
         human-picked (reduce (fn[m [k v]]
                             (assoc m k (inc (get m k 0))))
                           {} candidates)
         picked (count human-picked)

         good  (if good
                 good
                 (-> (get-good-candidates nm-re-coll cutpoint :ends ends)
                     first))
         good-size (count good)

         true+ (count (reduce (fn[m [k v]]
                                (if (get human-picked k)
                                  (assoc m k v) m))
                              {} good))
         false+ (- good-size true+)
         false- (- picked true+)]
     [:total total
      :cutpoint cutpoint
      :goodsize good-size
      :true+ true+
      :true+picked (float (/ true+ picked))
      :goodtrue+ (float (/ true+ good-size))
      :false+ false+
      :false- false-
      :false-picked (float (/ false- picked))
      :picked picked
      :dups (filter (fn[[k v]] (> v 1)) human-picked)]))


(defn plot-hit-dists
  [prot-name cnt-seqs ctx-delta res cutpt nms-stods & {:keys [file]}]
  (let [xs (range cnt-seqs)
        out (if file (fs/join file (str prot-name "-hitonly.png")) :display)]
    (render-chart
     out
     xs (map second nms-stods)
     "Sequence" "Sq RE to Hybrid"
     (str prot-name " Hits only to Hybrid."
          " Res: " res ", Cutpt: " cutpt)
     :series-label "sq/hbrid-distance"
     :legend true)))

(defn plot-cand-dists
  [prot-name cnt-seqs ctx-delta res cutpt nms-stods & {:keys [file]}]
  (let [xs (range cnt-seqs)
        out (if file (fs/join file (str prot-name "-candidates.png")) :display)]
    (render-chart
     out
     xs (map second nms-stods)
     "Sequence" "Sq RE to Hybrid"
     (str prot-name " Candidates to Hybrid."
          " Ctx Sz: " ctx-delta
          " Res: " res ", Cutpt: " cutpt)
     :series-label "sq/hbrid-distance"
     :legend true)))


(defn compute-candidate-info
  ""
  [sto-file candidate-file delta run
   & {:keys [refn xlate alpha crecut limit wz order
             plot-cre plot-dists]
      :or {refn DX||Y alpha (alphabet :rna) crecut 0.10 limit 15
           order :up}}]

  (let [[wz cres] (if wz
                    [wz]
                    (res-word-size sto-file :delta delta
                                   :xlate xlate :alpha alpha
                                   :limit limit :crecut crecut :cnt 5))
        [nm-re-sq pnm sz] (sto-re-dists wz candidate-file sto-file
                                        :refn refn :xlate xlate
                                        :delta delta :order order)
        ;; For hits, always use 0.9 on CDF
        Dy (if (zero? delta) 0.4 (case run 1 0.1, 2 0.0, -0.1))
        cutpt (select-cutpoint nm-re-sq :Dy Dy)
        [good bad] (get-good-candidates nm-re-sq cutpt)]
    (when (and plot-cre cres) (plot-cres cres))
    (when plot-dists
      (if (zero? delta)
        (plot-hit-dists pnm sz delta wz cutpt nm-re-sq :file plot-dists)
        (plot-cand-dists pnm sz delta wz cutpt nm-re-sq :file plot-dists)))
    [good bad cutpt nm-re-sq]))


(declare get-hitonly-final-names)

(defn compute-candidate-sets
  ""
  [sto-file candidate-file run delta
   & {:keys [refn xlate alpha crecut limit order
             cmp-ents plot-cre plot-dists]
      :or {refn DX||Y alpha (alphabet :rna) crecut 0.10
           limit 15 order :up}}]

  (let [[first-phase-entry-file final-entry-file]
        (get-hitonly-final-names candidate-file)

        final-entry-file (fs/fullpath final-entry-file)
        final-bad-entry-file (fs/replace-type
                              final-entry-file
                              (str "-bad." (fs/ftype final-entry-file)))
        [rna-only-good
         rna-only-bad
         cutpt1
         rna-nm-re-sq] (compute-candidate-info
                        sto-file candidate-file 0 run
                        :refn refn :xlate +RY-XLATE+ :alpha alpha
                        :limit limit :wz 6 :order order
                        :plot-cre plot-cre :plot-dists plot-dists)
        rna-only-perf-stats (when cmp-ents
                              (selection-perf
                               rna-nm-re-sq cmp-ents cutpt1
                               :good rna-only-good))
        _ (->> rna-only-good (map first)
               (#(gen-entry-file % first-phase-entry-file)))

        [good bad
         cutpt2
         nm-re-sq] (compute-candidate-info
                    sto-file first-phase-entry-file delta run
                    :refn refn :xlate xlate :alpha alpha
                    :crecut crecut :limit limit :order order
                    :plot-cre plot-cre :plot-dists plot-dists)
        perf-stats (when cmp-ents
                     (selection-perf
                      nm-re-sq cmp-ents cutpt2 :good good))]

    (->> good (#(gen-entry-nv-file % final-entry-file)))
    (->> bad (#(gen-entry-nv-file % final-bad-entry-file)))
    (entry-file->fasta-file final-entry-file)

    [[good bad cutpt2]
     [rna-only-good rna-only-bad cutpt1]
     rna-nm-re-sq
     nm-re-sq]))


(defn get-hitonly-final-names
  "Generate the canonical names for the output first phase 'hitonly',
   and 'final' entry files resulting from an FFP
   compute-candidate-sets selection from cmsearch targets.  CF is the
   candidate file name of the filtered results of the base cmsearch
   output.  Typically the csv file matching a cmsearch.out.
  "
  [cf]
  (let [suffix (first (re-find #"(\.[a-z]+|)\.cmsearch\.(csv|out|ent)$" cf))
        suffix-re (re-pattern suffix)
        hitonly (str/replace-re suffix-re "-hitonly.ent" cf)
        final (str/replace-re suffix-re "-final.ent" cf)]
    [hitonly final]))


(defn get-saved-ctx-size
  "If we have previously calculated the context size for sequences in
   STO, use that directly.  Previously calculated context sizes are
   recorded in sto on a #=GF line with subtype CTXSZ, followed by the
   size (an integer)
  "
  [sto]
  (let [l (->> sto io/read-lines
               (drop-until #(re-find #"^(#=GF\s+CTXSZ\s+[0-9]+|NC)" %))
               first (re-find #"CTXSZ\s+([0-9]+|)") second)]
    (if l (Integer. l) nil)))

(defn hit-context-delta
  [sto & {:keys [plot mindelta] :or {mindelta 200}}]
  (let [gfcsz (get-saved-ctx-size sto)]
    (if gfcsz
      [gfcsz true]
      (let [pts (xfold (fn[i]
                         (->> (compute-candidate-info
                               sto sto
                               (+ mindelta (* i 20)) 1
                               :refn jensen-shannon
                               ;;:xlate +RY-XLATE+ :alpha ["R" "Y"]
                               :crecut 0.01 :limit 19
                               :plot-dists false)
                              first (map second) mean))
                       (range 81))
            ms (map #(min %1 %2) pts (drop 1 pts))
            out (if plot
                  (->> sto fs/basename
                       (#(fs/replace-type % "-ctxsz.png"))
                       (fs/join plot))
                  :display)]
        (when plot
          (render-chart
           out
           (map #(+ mindelta (* 20 %)) (range (count ms))) ms
           "Size X 20" "RE/JSD" "RE to subseq"
           :series-label "Sub Seq Size"
           :legend true))
        [(+ mindelta (* 20 (first (pos (apply min pts) pts)))) false]))))

(defn hit-context-delta-db
  ""
  [sto & {:keys [gene cnt margin] :or {cnt 5 margin 0}}]
  {:pre [(string? gene)]}

  (let [sample (random-subset (read-seqs sto :info :name) cnt)
        nms-loci (map #(let [[nm [s e] sd] (entry-parts %)] [nm s e sd]) sample)
        features (map (fn[[nm s e sd]]
                        [[nm s e sd]
                         (edu.bc.bio.gaisr.db-actions/hit-features-query
                          nm s e)])
                      nms-loci)
        nms-lens
        (->> features
             (map (fn[[[nm s e sd :as entry] ctx]]
                    [entry
                     (ffirst
                      (keep (fn[m]
                              (let [nvs (m :nvs)
                                    x (keep #(when (and (= (% :name) "gene")
                                                        (= (% :value) gene))
                                               ;; -1 sd => hit end - ctx start
                                               ;;  1 sd => ctx end - hit start
                                               (let [{cs :start
                                                      ce :end
                                                      csd :strand}
                                                     (first (m :locs))
                                                     len (inc (if (= 1 csd)
                                                                (- ce s)
                                                                (- e cs)))]
                                               [len %]))
                                            nvs)]
                                (when (and (in (m :sftype) ["CDS" "gene"])
                                           (seq x))
                                  (first x))))
                            ctx))]))
             (filter second))
        base (+ margin
                (if (not (seq nms-lens))
                  (hit-context-delta sto :gene gene :cnt cnt)
                  (reduce (fn[sz [nm s]]
                            (first (div (+ sz s) 2)))
                          (-> nms-lens first second) (rest nms-lens))))]
    (* (floor (+ (/ base 100) 1)) 100)))




(defn krnn-seqs-clust
  ""
  [seqents & {:keys [xlate alpha delta crecut limit kinfo]
              :or {alpha (alphabet :rna)}}]
  {:pre [(or (and (vector? kinfo)
                  (reduce (fn[b x] (and b (integer? x))) true kinfo))
             (integer? kinfo)
             (nil? kinfo))]}

  (let [seqents (if (coll? seqents) seqents (read-seqs seqents :info :name))
        seqcnt (count seqents)
        kinfo (if (not (nil? kinfo)) kinfo [4 (long (->> seqcnt log2 ceil))])
        [kmin kmax] (if (vector? kinfo) kinfo [(min 4 (int (/ kinfo 2))) kinfo])

        ;; Get entries and their delta adjusted sequences and Goedel
        ;; number them for efficient map key access.
        entries  (map (fn[i es] [i es]) (iterate inc 0) seqents)
        seqs (get-adjusted-seqs (map second entries) delta)
        entries (into {} entries)
        coll (->> seqs
                  (map (fn[[e s]] [e (if xlate (seqXlate s :xmap xlate) s)]))
                  (mapv (fn[i es] [i es]) (iterate inc 0)))

        ;; Get the resolution window size for sequence words
        ;; (features), which controls vocabulary distribution content
        ;; (dictionaries), pmfs, and RE.
        [wz] (res-word-size (map first seqs) :delta delta
                            :xlate xlate :alpha alpha
                            :crecut crecut :limit limit)

        ;; Precompute all the resulting dictionaries of sequences with
        ;; word size.  This is a vector which is naturally indexed by
        ;; the numbering of the entries.
        coll-ffps (mapv (fn[[i [e sq]]] (probs wz sq)) coll)

        ;; Distance function that uses above precomputed distributions
        ;; as inputs keyed off entry/seq numbering.  Just call
        ;; coll-ffps on the numberings (i and j) of the inputs.  Each
        ;; arg has shape [i [ent-i sq-i]]; we only need i.
        distfn (fn[[x _x] [y _y]]
                 (jensen-shannon (coll-ffps x) (coll-ffps y)))

        ;; Compute distance matrix for knn graph computation.  keyfn
        ;; is first of input elements which gives the numbering for an
        ;; entry.  So, matrix is a map with keys [x y], x & y
        ;; associated numbers of coll entries
        keyfn first
        dm (clu/dist-matrix distfn coll :keyfn keyfn)


        ;; Distance function for S-Dbw cluster validity measure.
        ;; Every 'point' here should be word dictionary as seq-clus
        ;; (see below) is precomputed to be such.
        distfn2 (fn[l r] {:pre [(and (map? l) (map? r))]}
                  (jensen-shannon l r))

        ;; Averaging ('mean') function for S-Dbw validity measure.
        ;; This computes the 'hybrid' (minimized entropy) dictionary -
        ;; a 'centroid' dictionary.
        avgfn (fn
                ([sqs]
                   (hybrid-dictionary wz sqs))
                ([x & xs]
                   (hybrid-dictionary wz (cons x xs))))

        kcoll (map keyfn coll)]

    (println :wz wz"\n":kmax kmax)
    (for [k (range kmin kmax)
          :let [[krnngrph rnncntM knngrph]
                (clu/krnn-graph k #(get dm [%1 %2]) kcoll)]]
      (let [_ (println :start-clustering :k k (str-date))
            clusters (->> (clu/split-krnn k krnngrph rnncntM knngrph)
                          ;;(#(do (prn :G>k/G<k %) %))
                          (map #(gr/tarjan (keys %) %))
                          ;;(#(do (prn :sccs %) %))
                          (apply clu/refoldin-outliers krnngrph)
                          vec)
            ent-clus (let [coll (into {} coll)]
                       (map (fn[scc]
                              (map (fn[x]
                                     [(entries x) (-> x coll first)])
                                   scc))
                            clusters))
            seq-clus (let [coll (into {} coll)]
                       (map (fn[scc]
                              (let [v (xfold
                                       (fn[x] (coll-ffps x))
                                       (vec scc))]
                                [(avgfn v) v]))
                            clusters))
            _ (println :end-clustering :start-S-Dbw (str-date))]
        [(clu/S-Dbw-index distfn2 (vec seq-clus) :avgfn avgfn)
         ent-clus k]))))

(defn split-sto
  ""
  [stofile & {:keys [delta xlate alpha crecut limit kinfo]
              :or {alpha (alphabet :rna)}}]
  (let [basedir (fs/dirname stofile)
        name (->> stofile fs/basename (str/split #"-") first)
        clu-dir (fs/join basedir (str "CLU-" name))
        _ (when (not (fs/exists? clu-dir)) (fs/mkdir clu-dir))
        clud-dir (fs/join clu-dir (str "d" delta))
        _ (when (not (fs/exists? clud-dir)) (fs/mkdir clud-dir))
        clu-info (->> (krnn-seqs-clust
                       stofile
                       :delta delta
                       :xlate xlate :alpha alpha
                       :crecut crecut :limit limit
                       :kinfo kinfo)
                      (sort-by first <))
        fs (->> clu-info
                first second (map #(map first %))
                (map (fn[i ents]
                       (gen-entry-file
                        ents (fs/join clud-dir (str "clu-" i ".ent"))))
                     (iterate inc 1))
                doall)]
    [(map (fn[[s ents k]] [s k (count ents) (map count ents)]) clu-info) fs]))



(comment


(def L20-auto
     (compute-candidate-sets
      "/home/kaila/Bio/STOfiles/031113/AUTO/B/L20-auto-3.sto"
      "/home/kaila/Bio/STOfiles/031113/AUTO/B/CSV/Pass-3/L20-auto-3.sto.Assortprot1.fna.cmsearch.csv"
      3 1200
      :refn jensen-shannon
      :xlate +RY-XLATE+ :alpha ["R" "Y"]
      :crecut 0.01 :limit 19
      :plot-dists true))



(def bigrun2-info
     (compute-candidate-sets
      "/home/kaila/Bio/Tests/JSA/RNA_00012.sto"
      "/home/kaila/Bio/Tests/RNA_00012.sto.Assortprots2.fna.cmsearch.csv"
      "/home/kaila/Bio/Tests/RNA_00012.sto.Assortprots2-hitonly.ent"
      "/home/kaila/Bio/Tests/RNA_00012.sto.Assortprots2-final.ent"
      (hit-context-delta "/home/kaila/Bio/Tests/RNA_00012.sto" :gene "rpsO")
      :refn jensen-shannon :xlate +RY-XLATE+ :alpha ["R" "Y"] :crecut 0.01 :limit 19
      :order :up :ends :low
      :plot-cre false :plot-dists true :cutoff 33/100))

(def rfam-00504-AP2-test
     (compute-candidate-sets
      "/home/jsa/Bio/FreqDicts/NewRFAM/STOS/RF00504-seed-NC.sto"
      "/home/jsa/Bio/FreqDicts/NewRFAM/STOS/CSV/RF00504-seed-NC.sto.Assortprots2.fna.cmsearch.csv"
      "/home/jsa/Bio/FreqDicts/NewRFAM/STOS/RF00504-AP2-hitonly.ent"
      "/home/jsa/Bio/FreqDicts/NewRFAM/STOS/RF00504-AP2-final.ent"
      (hit-context-delta "/home/jsa/Bio/FreqDicts/NewRFAM/STOS/RF00504-seed-NC.sto" :gene "gcvT")
      :refn jensen-shannon :xlate +RY-XLATE+ :alpha ["R" "Y"] :crecut 0.01 :limit 19
      :order :up :ends :low
      :plot-cre false :plot-dists true :cutoff 37/100))

)



#_(def sample-metag-cres
     (->> (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")
          first read-seqs
          (take 20)
          (map #(pxmap (fn[k] [k (CREl k %) (count %)]) 2 (range 3 16)))))

#_(doseq [cres (take 3 sample-metag-cres)]
  (let [ks (map first cres)
        cnt (third (first cres))
        cres (map second cres)]
    (incanter.core/view
     (incanter.charts/scatter-plot
      ks cres
      :x-label "feature length"
      :y-label "Commulative Relative Entropy"
      :title (str "CRE(l, F/F^), len: " cnt
                  " Fmax: " (round (ln cnt)))))))

#_(doseq [cres foocres]
  (let [ks (map first cres)
        cnt (third (first cres))
        cres (map second cres)]
    (incanter.core/view (incanter.charts/scatter-plot
                         ks cres
                         :x-label "feature length"
                         :y-label "Commulative Relative Entropy"
                         :title (str "CRE(l, F/F^), len: " cnt
                                     " Fmax: " (round (ln cnt)))))))

#_(def count-map
     (reduce (fn[m file]
               (reduce (fn[m sq]
                         (let [cnt (count sq)]
                           (assoc m cnt (inc (get m cnt 0)))))
                       m (read-seqs file)))
             {} (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")))


#_(def count-maps
       (map #(do [(fs/basename %)
                  (reduce (fn[m sq]
                            (let [cnt (count sq)]
                              (assoc m cnt (inc (get m cnt 0)))))
                          {} (read-seqs %))])
            (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")))

#_(sort-by
   first
   (map #(do [(first %) (->> % second (sort-by first >)
                             ((fn[sq]
                                [(take 5 sq) (drop (- (count sq) 5) sq)])))])
        count-maps))

