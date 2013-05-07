;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                               O P E R O N S                              ;;
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

(ns edu.bc.bio.rdb.operons

  "Build operon table and connections into std biosql tables.  This is
   typically a load once and use given a snapshot of the operon data."

  (:require [clojure.contrib.sql :as sql]
            [org.bituf.clj-dbcp :as dbcp]
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.contrib.io :as io]
            [clojure.zip :as zip]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        [edu.bc.log4clj :only [create-loggers log>]]
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)])

  (:import javax.sql.DataSource
           com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource))


(create-loggers
 [[:logger "rdbLogger"
   {:level :info :appenders "rdbinfo"}]
  [:appender "rdbinfo"
   {:type :rolling
    :filespec "./Logs/rdbinfo.log"
    :max-version 4
    :max-size "5MB"}]])


(def db-pw (atom ""))

;;; Define MySql datasource connection pool.  This affords reasonable
;;; (for mysql...) parallelization of concurrent queries and keeps the
;;; ansync web service reasonably happy (though mysql is the bottleneck...
;;;
(def mysql-ds
     (dbcp/db-spec
      (let [ds (dbcp/mysql-datasource
                "127.0.0.1:3306" "biosql" "root" @db-pw)]
        (dbcp/set-max-active! ds 50)
        (dbcp/set-min-max-idle! ds 5 20)
        ds)))



(defn get-line-info []
  (let [l (read-line)]
    (when (not (nil? l))
      (let [parts (vec (str/split #"\t" l))]
        (if (< (count parts) 10)
          (do (log> "rdbLogger" :warn
                    "Missing Operon Info: '~A' --> '~A'" l parts)
              (get-line-info))
          (let [name (parts 0)
                opid (Integer. (parts 1))
                loc (parts 2)
                loc (str/split #"\.\." loc)
                start (Integer. (first loc))
                end (Integer. (second loc))
                strand (parts 3)
                strand (if (= strand "+") 1 -1)
                pid (parts 5)
                cog (parts 9)
                opr-info [opid start end strand pid cog]]
            [name opid opr-info]))))))


(defn bioentry-id [gbname]
  (let [fields "bioentry_id"
        table "bioentry as be"
        constraint (str "be.name=" \" gbname \")
        stmt (str "select " fields
                  " from " table
                  " where " constraint)
        result (first (sql/with-connection mysql-ds
                        (sql/with-query-results qresults [stmt]
                          (doall qresults))))]
    (result :bioentry_id)))


(defn print-update-operon-tables [gbname opr-info]
  (when (not (empty? opr-info))
    (let [beid (bioentry-id gbname)
          opids (sort (reduce #(conj %1 (first %2)) #{} opr-info))
          opr-rows (map #(do %1 [beid]) opids)
          opr-loc-rows (sort-by first opr-info)]
      (println "***" gbname)
      (println opr-rows)
      (println opr-loc-rows))))


(defn update-operon-tables [gbname opr-info]
  (log> "rdbLogger" :info
        "~A operons: ~A" gbname (count opr-info))
  (println gbname "-->" (count opr-info))
  (when (not (empty? opr-info))
    (let [beid (bioentry-id gbname)
          opids (sort (reduce #(conj %1 (first %2)) #{} opr-info))
          opr-rows (repeat (count opids) [beid])
          opr-loc-rows (sort-by first opr-info)]
      (sql/with-connection mysql-ds
        (sql/transaction
         (apply sql/insert-values
                :operon [:bioentry_id] opr-rows)
         (apply sql/insert-values
                :operon_loc
                [:operon_id :start_pos :end_pos :strand :protein_id :cog]
                opr-loc-rows))))))


;;;"/data2/BioData/test.txt"
(defn populate-operon-tables [file]
  (with-open [opr (io/reader file)]
    (binding [*in* opr]
      (read-line) ; toss header
      (with-local-vars [info [], gbname "", opr-loc-set #{}]
        (while  (let [x (get-line-info)] (var-set info x))
          (let [[nm opid opinfo] (var-get info)
                gbnm (var-get gbname)]
            (if (not= gbnm nm)
              (do (update-operon-tables gbnm (var-get opr-loc-set))
                  (var-set opr-loc-set #{opinfo})
                  (var-set gbname nm))
              (let [op-set (var-get opr-loc-set)
                    op-set (conj op-set opinfo)]
                (var-set opr-loc-set op-set)))))
        (update-operon-tables (var-get gbname) (var-get opr-loc-set))))))




(def counter (atom 0))

(defn log-info [gbname opr-info]
  (when (not (empty? opr-info))
    (let [opids (sort (reduce #(conj %1 (first %2)) #{} opr-info))]
      (swap! counter #(+ % (count opids)))
      (when (not= @counter (last opids))
        (raise :type :bad-operon-id :opids opids :count @counter))
      (println gbname "-->" (count opr-info)))))

(defn check-operon-file [file]
  (swap! counter #(do % 0))
  (with-open [opr (io/reader file)]
    (binding [*in* opr]
      (read-line) ; toss header
      (with-local-vars [info [], gbname "", opr-loc-set #{}]
        (while  (let [x (get-line-info)] (var-set info x))
          (let [[nm opid opinfo] (var-get info)
                gbnm (var-get gbname)]
            (if (not= gbnm nm)
              (do (log-info gbnm (var-get opr-loc-set))
                  (var-set opr-loc-set #{opinfo})
                  (var-set gbname nm))
              (let [op-set (var-get opr-loc-set)
                    op-set (conj op-set opinfo)]
                (var-set opr-loc-set op-set)))))
        (log-info (var-get gbname) (var-get opr-loc-set))))))
