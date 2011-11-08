;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                             N E T - U T I L S                            ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011 Trustees of Boston College                            ;;
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

(ns edu.bc.net-utils

  "General Network utility functions and macros.  Basically these
   resources are fairly general and intended to be usable on most any
   project that needs to talk http, ftp, etc."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.contrib.properties :as prop]
            [clojure.xml :as xml]
            [clj-http.client :as client]
            [net.cgrand.enlive-html :as html]

            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)])

  (:import (org.apache.commons.net.ftp
            FTP FTPClient FTPFile FTPListParseEngine)
           (java.net URL)))


;;; (defn refseq-filter []
;;;   (reify
;;;    FTPFileFilter
;;;    (^boolean accept [this ^FTPFile ftpfile]
;;;          (println "ACCEPT was called with: " (.getName ftpfile))
;;;          (if (not (.isFile ftpfile))
;;;            false
;;;            (let [name (.getName ftpfile)]
;;;              (if (re-find #"[0-9]+\.genomic\.(gbff|fna)\.gz" name)
;;;                true
;;;                false))))))
;;;
;;; (def *rsf* (refseq-filter))
;;; (instance? FTPFileFilter *rsf*)




(defmulti open "Open an ftp connection." class)


(defmethod open String [s]
  (open (java.net.URL. s)))


(defmethod open URL [url]
  (let [client (FTPClient.)]
    (.connect client
              (.getHost url)
              (if (= -1 (.getPort url))
                (.getDefaultPort url)
                (.getPort url)))
    client))




(defmacro with-ftp [[client url & extra-bindings] & body]
  `(let [u# (URL. ~url)
         ~client (open u#)
         res# (atom nil)
         ~@extra-bindings]
     (try
       (do
         (if (.getUserInfo u#)
           (let [[uname# pass#] (.split (.getUserInfo u#) ":" 2)]
             (.login ~client uname# pass#)))
         (assert (.changeWorkingDirectory ~client (.getPath u#)))
         (.setFileType ~client FTP/BINARY_FILE_TYPE)
         (reset! res# (do ~@body))
         @res#)
       (finally
        (.disconnect ~client)))))


(defmacro with-paged-ftp [[client url & extra-bindings] & body]
  `(let [u# (URL. ~url)
         ~client (open u#)
         res# (atom [])
         page-size# 50 ; Need to support this in extra-bindings!
         ~@extra-bindings]
     (try
       (do
         (if (.getUserInfo u#)
           (let [[uname# pass#] (.split (.getUserInfo u#) ":" 2)]
             (.login ~client uname# pass#)))
         (assert (.changeWorkingDirectory ~client (.getPath u#)))
         (.setFileType ~client FTP/BINARY_FILE_TYPE)
         (let [pwd# (.printWorkingDirectory ~client)
               engine# (.initiateListParsing ~client pwd#)]
           (while (. engine# hasNext)
             (let [~'$chunk (. engine# getNext page-size#)]
               (reset! res# (concat @res# (do ~@body))))))
         @res#)
       (finally
        (.disconnect ~client)))))




(defn only-dirs [ftp-set]
  (filter #(.isDirectory %) ftp-set))


(defn only-files [ftp-set]
  (filter #(.isFile %) ftp-set))


(defn file-names [ftp-set]
  (keep #(when (.isFile %) (.getName %)) ftp-set))




(defn ftp-dir
  ([url]
     (with-ftp [client url]
       (.listFiles client)))
  ([url fn]
     (with-paged-ftp [client url]
       (fn $chunk))))


(defn list-files [url]
  (ftp-dir url file-names))


(defn list-dirs [url]
  (ftp-dir url (fn [ftp-set] (map #(.getName %) (only-dirs ftp-set)))))


(defn retrieve-file [url fname & [local-file]]
  (with-ftp [client url]
    (if local-file
      (with-open [outstream (java.io.FileOutputStream. local-file)]
        (.retrieveFile client fname outstream))
      (.retrieveFile client fname))))


;;; --------------- HTTP Client ------------

(defn fetch-url [url]
  (html/html-resource (java.net.URL. url)))

(defn fetch-elements-by-tag [url-map tags]
  (map html/text (html/select url-map tags)))

(defn http-binary-download [url local-file]
  (with-open [is (io/input-stream url)
              os (io/output-stream local-file)]
    (io/copy is os)))

;(:body (client/get "http://www.bork.embl.de/~arumugam/Qin_et_al_2010/"))

