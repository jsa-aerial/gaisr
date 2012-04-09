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


(defn parse-config-file [filespec]
  (let [lseq (filter #(not (or (= "" %) (re-find #"^#+" %)))
		     (io/read-lines filespec))]
    (reduce
     (fn[m l]
       (condp #(re-find %1 %2) l
	 #"^Hitfile-Dir *:"
	 (assoc m :hitfile-dir (second (str/split #": *" l)))
	 
	 #"^CM-Dir *:"
	 (assoc m :cmdir (second (str/split #": *" l)))
	 
	 #"^CM-Search +hitfile *:"
	 (assoc m :hitfile (second (str/split #": *" l)))
	 
	 #"^CM-Search *:"
	 (assoc m :hitfile nil)

	 (let [hitfile (m :hitfile)
	       [cm hf] (str/split #" *, *" l)
	       the-hitfile (or hitfile hf)]
	   (assert (not (and hitfile hf)))
	   (let [cms-map (get m :cmsearchs {})
		 cm-hfs (get cms-map cm [])
		 cm-hfs (conj cm-hfs the-hitfile)]
	     (assoc m :cmsearchs (assoc cms-map cm cm-hfs))))))
     {} lseq)))