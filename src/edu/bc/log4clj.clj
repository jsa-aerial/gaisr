;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                          L O G 4 C L J                                   ;;
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

(ns edu.bc.log4clj
  "Attempt to make some sort of reasonable abstraction over log4j.  Or
   maybe we will give up and just reimplement log4cl in clj.  Java
   logging a gigantic awful mess of over complexity to attain what is
   actually pretty simple and straightforward.  Log4cl does everything
   log4j does, while also being more extensible and dynamic, yet is
   about 1/10 the complexity.  Clojure, being JVM based, seems like it
   would be reasonable to just use log4j, but the underlying lib
   complexity is nuts."
  (:require [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.properties :as prop]
            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        [clojure.pprint
         :only [cl-format]])

  (:import (java.util Properties Vector)
           ;; Useless... (org.apache.commons.logging Log)
           (org.apache.log4j PropertyConfigurator
                             LogManager
                             Logger)))


;;; Logging levels for log4j and with c.c.logging.  As with
;;; c.c.logging these are case as keywords as that is saner
;;;
;;; lowest ... highest
;;;
;;; :ALL :TRACE :DEBUG :INFO :WARN :ERROR :FATAL :OFF


(defparameter *levels*
     {:all   org.apache.log4j.Level/ALL
      :trace org.apache.log4j.Level/TRACE
      :debug org.apache.log4j.Level/DEBUG
      :info  org.apache.log4j.Level/INFO
      :warn  org.apache.log4j.Level/WARN
      :error org.apache.log4j.Level/ERROR
      :fatal org.apache.log4j.Level/FATAL
      :off   org.apache.log4j.Level/OFF})

(defparameter *appenders*
     {:console "org.apache.log4j.ConsoleAppender"
      :file    "org.apache.log4j.FileAppender"
      :rolling "org.apache.log4j.RollingFileAppender"})

(defparameter *options*
     {:level       "Level" ; Don't need now, but will when reimpl
      :name        "Name"  ;     "
      :layout      "layout.ConversionPattern"
      :pat         "org.apache.log4j.PatternLayout"
      :filespec    "File"
      :buf-size    "BufferSize"
      :buf-io      "BufferedIO"
      :append      "Append"
      :max-size    "MaxFileSize"
      :max-version "MaxBackupIndex"})


(defn config-logging
  "Initial startup logging configuration.  Uses property file
   PROP-DESIGNATOR, or if not given, defaults to log4j.properties.
   ***NOTE: more or less obsolete now we have configure-logging."
  ([]
     (config-logging "log4j.properties"))
  ([prop-designator]
     (if (string? prop-designator)
       (PropertyConfigurator/configure
        (fs/fullpath prop-designator))
       (PropertyConfigurator/configure
        prop-designator))))


;;; Not really useful.  The Commons logging wrapper sucks anyway and
;;; since this is log4clj, it is already more or less wedded to either
;;; 1. log4j or 2. its own direct implementation of more or less what
;;; log4j is/does.  Currently it is 1.  Either way, the commons access
;;; api is confusing and a useless/worthless layer.
;;; (defn get-log
;;;   ([]
;;;      (get-log (str *ns*)))
;;;   ([name]
;;;      (LogFactory/getLog name)))

(defn get-logger
  "Return the logger named NAME (a string or a logger) or if NAME is not given,
   use the name of the current name space.  If NAME is a logger, just
   return it."
  ([]
     (get-logger (str *ns*)))
  ([name]
     (if (string? name)
       (if (= name "root")
         (LogManager/getRootLogger)
         (LogManager/getLogger name))
       name)))

(defn get-root-logger []
  (get-logger "root"))


(defn get-all-loggers
  "Return a vector of all the current loggers"
  []
  (loop [jloggers (LogManager/getCurrentLoggers)
         loggers [(LogManager/getRootLogger)]]
    (if (not (. jloggers hasMoreElements))
      loggers
      (recur jloggers
             (conj loggers (. jloggers nextElement))))))


(defn all-loggers-map
  "Return a map of all current loggers indexed by their name"
  []
  (reduce (fn [m l] (assoc m (. l getName) l))
          {} (get-all-loggers)))


(defn logger-appenders
  "Return all the appenders for logger LOGGER, a string naming the
  logger or a logger"
  [logger]
  (let [logger (get-logger logger)]
    (loop [jappenders (. logger getAllAppenders)
           appenders []]
      (if (not (. jappenders hasMoreElements))
        appenders
        (recur jappenders
               (conj appenders (. jappenders nextElement)))))))


(defn all-appenders
  "Return all the current appenders.  This set is the union of all the
   appenders of all the current loggers"
  []
  (apply set/union
         (map logger-appenders (get-all-loggers))))


(defn close-appender
  "Close the appender A (an appender as obtained from logger-appenders
   or all-appenders).  Once, closed, this appender instance is no
   longer useable, however it is always possible to reconfigure
   logging with corresponding new appenders - with same names and on
   same loggers"
  [a]
  (doto a .close))


(defn close-all-appenders
  "For all A in (all-appenders), (close-appender A)"
  []
  (doseq [a (all-appenders)]
    (close-appender a)))



;;; (logger-decl "refseq-ftp" {:level :info :appenders "netutils"})
;;; (logger-decl
;;;  "rootLogger"
;;;  {:level :info :appenders "console"})
;;;
(defn logger-decl [name attrs]
  (let [lname (if (not= name "rootLogger")
                (str "log4j.logger." name)
                (str "log4j." name))
        level (. (*levels* (get attrs :level :info)) toString)
        additivity (str "log4j.additivity." name)
        appenders (get attrs :appenders "console")
        appenders (if (coll? appenders)
                    (str/join ", " appenders)
                    appenders)]
    [lname (str/join ", " [level appenders])
     (list [additivity "false"])]))


(defn file-attrs [type attrs]
  (if (= type :console)
    []
    (let [buf-io (get attrs :buf-io "true")
          buf-size (get attrs :buf-size "1024")
          filespec (get attrs :filespec (str "./" (ns-name *ns*) ".log"))
          append (get attrs :append "true")]
      (map
       (fn [k v]
         [(str "." (get *options* k)) v])
       [:buf-io :buf-size :filespec :append]
       [buf-io buf-size filespec append]))))

(defn rolling-attrs [type attrs]
  (if (not= :rolling type)
    []
    (let [max-size (str (get attrs :max-size "5MB"))
          max-version (str (get attrs :max-version 2))]
      (map
       (fn [k v]
         [(str "." (get *options* k)) v])
       [:max-size :max-version]
       [max-size max-version]))))


(defn appender-decl [name attrs]
  (let [name (str "log4j.appender." name)
        config-type (get attrs :type :console)
        config-type (if (*appenders* config-type)
                      config-type
                      :console)
        type (*appenders* config-type)

        layout-key ".layout"
        layout-val (get *options* :pat)
        convpat-key ".layout.ConversionPattern"
        convpat-val (get attrs :pat "%-5p %c %d{DATE}: %m%n")
        base-attrs [[layout-key layout-val]
                    [convpat-key convpat-val]]

        file-attrs (file-attrs config-type attrs)
        rolling-attrs (rolling-attrs config-type attrs)]
    [name type (concat base-attrs file-attrs rolling-attrs)]))


(defn javaize-log-decl [decl props]
  (let [[kind name attrs] decl]
    (assert (#{:logger :appender} kind))
    (if (= kind :logger)
      (let [[name kind attrs] (logger-decl name attrs)]
        (.setProperty props name kind)
        (doseq [[k v] attrs]
          (.setProperty props k v)))
      ;; Need to refactor now that loggers have attrs...
      (let [[name kind attrs] (appender-decl name attrs)]
        (.setProperty props name kind)
        (doseq [[k v] attrs]
          (.setProperty props (str name k) v)))))
  props)


(defn make-log-props [decls]
  (let [log-props (Properties.)]
    (doseq [decl decls]
      (javaize-log-decl decl log-props))
    log-props))


(defn configure-logging [decls]
  (PropertyConfigurator/configure
   (make-log-props decls)))


(defn activate-logging []
  (doseq [a (all-appenders)]
    (. a activateOptions)))


(defn create-loggers [decls]
  (configure-logging decls)
  (activate-logging))



;;; (create-loggers
;;;  [[:logger "rootLogger"
;;;    {:level :info :appenders "console"}]

;;;   [:appender "console"
;;;    {:type :console :pat "%-5p %d{DATE}: %m%n"}]

;;;   [:logger "refseq-ftp"
;;;    {:level :info :appenders "netutils"}]

;;;   [:appender "netutils"
;;;    {:type :rolling
;;;     :filespec "./Logs/netutils.log"
;;;     :max-version 4
;;;     :max-size 1000 ; TESTING!!
;;;     }]
;;;   ])




(defn get-level
  "Return the logging level for LOGGER.  If not given, the default is
   to try a logger for the current name space name and then a parent
   for that.  If LOGGER is a logger just use it, if it is a string
   use (get-logger LOGGER)"
  ([]
     (keyword
      (..
       (or (.. (get-logger) getLevel)
           (.. (get-logger) getParent getLevel))
       toString toLowerCase)))
  ([logger]
     (let [logger (get-logger logger)]
       (keyword (.. logger getLevel toString toLowerCase)))))


(defn set-level!
  "Sets the logger level.  If not given, logger defaults to the root
   logger.  Levels are given as keywords and are from lowest to
   highest: :all :trace :debug :info :warn :error :fatal :off"
  ([level]
     (let [root-logger (Logger/getRootLogger)]
       (set-level! root-logger (*levels* level))))
  ([logger level]
     (.setLevel logger (*levels* level))))


(defn set-debug!
  "Sets log level of LOGGER to :debug.  If LOGGER is not given
  defaults to the root logger"
  ([]
     (set-level! :debug))
  ([logger]
     (set-level! logger :debug)))

(defn set-info!
  "Sets log level of LOGGER to :info.  If LOGGER is not given
  defaults to the root logger"
  ([]
     (set-level! :info))
  ([logger]
     (set-level! logger :info)))




(defn log>
  "Request LOGGER to log a message of level LEVEL, with message given
   as the resulting string of applying format string to ARGS (as
   defined by cl-format"
  [logger level message & args]
  (assert (*levels* level))
  (let [logger (get-logger logger)]
    (. logger log
       (*levels* level)
       (apply cl-format nil message args))))

