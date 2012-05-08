;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                             J O B - C O N F I G                          ;;
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

(ns edu.bc.bio.job-config

  "Gaisr pipeline job configuration."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clj-shell.shell :as sh]
            [edu.bc.fs :as fs])
  (:use clojure.contrib.math
        edu.bc.utils
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        ))


(defn process-sto-files [m l]
  (let [sto (str/trim l)]
    (let [stos (conj (get m :cmbuilds []) sto)]
      (assoc m :cmbuilds stos))))

(defn process-cm-files [m l]
  (let [sto (str/trim l)]
    (let [stos (conj (get m :calibrates []) sto)]
      (assoc m :calibrates stos))))

(defn process-hit-files [m l]
  (let [hfg (m :hitfile)
        [cm hfs] (str/split #" *, *" l)
        hitfile (or hfg hfs)]
    (assert (not (and hfg hfs)))
    (let [xref-map (get m :cmsearchs {})
          hf-cms (get xref-map hitfile [])
          hf-cms (conj hf-cms cm)]
      (assoc m :cmsearchs (assoc xref-map hitfile hf-cms)))))


(defn parse-config-file [filespec]
  (let [lseq (filter #(not (or (= "" (str/trim %))
                               (re-find #"^ *#+" %)))
                     (io/read-lines filespec))
        add-comment (fn [m l]
                      (assoc m :comments
                             (conj (get m :comments []) l)))
        directive (fn [m d]
                    (assoc m :directive d))
        getf (fn[m k l]
               (assoc m k (str/trim
                           (second (str/split #": *" l)))))]
    (reduce
     (fn[m l]
       (condp #(re-find %1 %2) l
         #"^ *#+" ; NOTE: need to remove re-find from lseq filter
         (add-comment m l)

         #"^Hitfile-Dir *:"
         (getf m :hitfile-dir l)

         #"^CM-Dir *:"
         (getf m :cmdir l)

         #"^STO-Dir *:"
         (getf m :stodir l)

         #"^CM-Build[ :]*"
         (directive (assoc m :cmbuild (m :cmdir)) :cmbuild)

         #"^CM-Calibrate[ :]*"
         (directive (assoc m :cmcalibrate (m :cmdir)) :calibrate)

         #"^CM-Search +hitfile *:"
         (directive
          (assoc m :hitfile (second (str/split #": *" l)))
          :cmsearch)

         #"^CM-Search *:"
         (directive (assoc m :hitfile nil) :cmsearch)

         #"^Gen-CSVs[ :]*"
         (directive
          (assoc m :gen-csvs
                 (get (getf m :csvdir l)
                      :csvdir
                      (m :cmdir)))
          :gen-csvs)

         (case (m :directive)
               :cmbuild (process-sto-files m l)
               :calibrate (process-cm-files m l)
               :cmsearch (process-hit-files m l)
               :gen-csvs m)))
     {} lseq)))

