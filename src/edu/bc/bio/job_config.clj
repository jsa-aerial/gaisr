;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                             J O B - C O N F I G                          ;;
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

(ns edu.bc.bio.job-config

  "Gaisr pipeline job configuration."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        ))

;;; Near term extensions:
;;;
;;;  1. Regexs for file denoting.  So, instead of listing all cms to
;;;     run against a hitfna, for example, you could give *.cm for a
;;;     CM-Search directive.  Extend other action directives (cmbuild,
;;;     etc.) similarly.
;;;
;;;  2. Include options for action directives.  In particular, EV for
;;;  CM-Search and RW (rewrite


(defn process-1-file [m directive file]
  (let [files (conj (get m directive []) file)]
    (assoc m directive files)))


(defn process-sto-files [m l]
  (let [stodir (m :stodir)]
    (reduce (fn[m sto] (process-1-file m :cmbuilds sto))
            m (fs/glob (fs/fullpath (fs/join stodir l))))))

(defn process-cm-files [m l]
  (process-1-file m :calibrates l))


(defn process-hit-files [m l]
  (let [hfg (m :hitfile)
        [cm hfs] (str/split #"\s*,\s*" l)
        hitfile (or hfg hfs)]
    (assert (not (and hfg hfs)))
    (let [hitdir (m :hitfile-dir)
          hitfs (fs/re-directory-files hitdir (or hfg hfs))]
      (reduce (fn[m hf]
                (let [xref-map (get m :cmsearchs {})
                      hf-cms (get xref-map hf [])
                      hf-cms (conj hf-cms cm)]
                  (assoc m :cmsearchs (assoc xref-map hf hf-cms))))
              m hitfs))))

(defn process-fna-files [m l]
  (reduce (fn[m fna] (process-1-file m :blast fna))
          m (fs/glob (fs/fullpath l))))

(defn process-sccs-subline [m l]
  (let [[d v] (str/split #"\s*(:| )\s*" l)]
    (assoc m :sccs (conj (get m :sccs []) [(keyword d) v]))))

(defn process-gen-csv-subline [m l]
  (let [[d v] (str/split #"\s*(:| )\s*" l)
        d (keyword d)
        v (if (and (= d :aggregate) (= v "csv-dir")) (first (m :gen-csvs)) v)]
    (assoc m :gen-csvs (conj (get m :gen-csvs []) [(keyword d) v]))))

(defn process-gen-sto-subline [m l]
  (let [[d v] (str/split #"\s*(:| )\s*" l)]
    (assoc m :gen-stos (conj (get m :gen-stos []) [(keyword d) v]))))

(defn process-directive-options [m d l]
  (let [args (->> l (str/split #"\s*:\s*") second)
        args (if (not args)
               []
               (->> args
                    (str/split #"\s*,\s*")
                    (map #(->> % (str/split #"\s+")
                               ((fn[[k v]]
                                  (let [v (if (= k "dir")
                                            v
                                            (->> v read-string
                                                 ((fn[x] (if (symbol? x)
                                                           (keyword (name x))
                                                           x)))))]
                                    [(keyword k) v])))))
                    flatten
                    vec))]
    (assoc m d (conj (get m d []) [:args args]))))

(defn parse-config-file [filespec]
  (let [lseq (filter #(not (or (= "" (str/trim %))
                               (re-find #"^\s*#+" %)))
                     (io/read-lines filespec))
        add-comment (fn [m l]
                      (assoc m :comments
                             (conj (get m :comments []) l)))
        directive (fn [m d]
                    (assoc m :directive d))
        getf (fn[m k l]
               (assoc m k (str/trim
                           (second (str/split #":\s*" l)))))]

    ;; Basically a simple state machine.  Not quite an FSA as there is
    ;; some memory involved. But no stack either.
    (reduce
     (fn[m l]
       (condp #(re-find %1 %2) l
         #"^\s*#+" ; NOTE: need to remove re-find from lseq filter
         (add-comment m l)

         #"^BLAST\s*:"
         (directive (process-directive-options m :blast l) :blast)

         #"^Hitfile-Dir\s*:"
         (getf m :hitfile-dir l)

         #"^CM-Dir\s*:"
         (getf m :cmdir l)

         #"^STO-Dir\s*:"
         (getf m :stodir l)

         #"^CM-Build[\s:]*"
         (directive (assoc m :cmbuild (m :cmdir)) :cmbuild)

         #"^CM-Calibrate[\s:]*"
         (directive (assoc m :cmcalibrate (m :cmdir)) :calibrate)


         ;; NEEDS TO BE FIXED: This clause must come before other
         ;; cmsearch clauses to ensure both eval and hitfile inclusion.
         #"^CM-Search\s+eval\s*:\s*[0-9]+[.0-9]*,\s*hitfile\s*:"
         (let [x (str/split #"," l)
               [ev hf] (map #(second (str/split #":\s*" %)) x)]
           (directive
            (assoc (assoc m :eval (Double. ev)) :hitfile hf)
            :cmsearch))

         #"^CM-Search\s+hitfile\s*:"
         (directive
          (assoc m :hitfile (second (str/split #":\s*" l)))
          :cmsearch)

         #"^CM-Search\s+eval\s*:"
         (directive
          (assoc m :eval (Double. (second (str/split #":\s*" l))))
          :cmsearch)

         #"^CM-Search\s*:"
         (directive (assoc m :hitfile nil) :cmsearch)


         #"^Gen-CSVs\s+spread\s*:*"
         (let [x (str/split #"," l)
               [sp dir] (map #(str/trim (second (str/split #":\s*" %))) x)]
           (directive
            (assoc (assoc m :spread (Integer. sp)) :gen-csvs dir)
            :gen-csv))

         #"^Gen-CSVs[\s:]*"
         ;; HACK, this needs to be folded into process-directive-options
         (directive
          (assoc m :gen-csvs [(get (getf m :csvdir l) :csvdir (m :cmdir))])
          :gen-csvs)

         #"^SCCS\s*:"
         (directive (process-directive-options m :sccs l) :sccs)

         #"Gen-Stos\s*:"
         (directive m :gen-stos)

         (case (m :directive)
               :blast (process-fna-files m l)
               :cmbuild (process-sto-files m l)
               :calibrate (process-cm-files m l)
               :cmsearch (process-hit-files m l)
               :gen-csvs (process-gen-csv-subline m l)
               :sccs (process-sccs-subline m l)
               :gen-stos (process-gen-sto-subline m l))))
     {} lseq)))

