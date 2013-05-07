;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                              T A X O N O M Y                             ;;
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

(ns edu.bc.bio.rdb.taxonomy
  (:require [clojure.contrib.sql :as sql]
            [org.bituf.clj-dbcp :as dbcp]
            [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        [edu.bc.log4clj :only [create-loggers log>]]
        [clojure.pprint
         :only [cl-format]])

  (:import javax.sql.DataSource
           com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource))

;;;(use 'edu.bc.utils :reload)


(def nodes (ref (hash-map)))

(defn process-taxonomy-nodes [file]
  (do-text-file [file]
    (let [[id pid] (str/split #"\t\|\t" $line)
          id  (with-in-str id (read))
          pid (with-in-str pid (read))]
      (dosync (commute nodes (fn [hm k v] (assoc hm k v)) id pid)))))


;;;(process-taxonomy-nodes "/data2/BioData/taxdata/nodes.dmp")


(defn ancestors [id]
  (loop [x id
         v (list)]
      (if (= x 1)
        v
        (let [p (@nodes x)]
          (recur p (conj v x))))))

;;;(ancestors 1279)

(defn create-node-ancestors-file [file]
  (io/with-out-writer (io/file-str file)
    (doseq [id (keys @nodes)]
      (println (str id ": " (str/join ", " (ancestors id)))))))




(def names (ref (hash-map)))

(defn process-taxonomy-names [file]
  (do-text-file [file]
    (let [[id name-txt _ name-class] (str/split #"\t\|\t" $line)
          name-class (first (str/split #"\t\|" name-class))]
      (when (= name-class "scientific name")
        (let [id  (with-in-str id (read))]
          (dosync (commute names (fn [hm k v] (assoc hm k v)) id name-txt)))))))

;;;(process-taxonomy-names "/data2/BioData/taxdata/names.dmp")


(defn ancestors-as-names [id]
  (loop [x id
         v (list)]
    (if (= x 1)
      v
      (let [p (@nodes x)
            n (@names x)]
        (recur p (conj v n))))))

;;;(ancestors-as-names 2)

(defn create-node-ancestor-names-file [file]
  (io/with-out-writer (io/file-str file)
    (doseq [id (keys @nodes)]
      (when (not= id 1) ; don't include 'root' node
        (println (str id "\t|\t" (str/join ", " (ancestors-as-names id))))))))

;;;(create-node-ancestor-names-file "/data2/BioData/taxdata/nodes-ancestor-names.txt")




;;; ------------------------------------------------------------------------;;;
;;; Put the stuff into the DB ....


(def db-pw (atom ""))

;;; Define MySql datasource connection pool.  This affords reasonable
;;; (for mysql...) parallelization of concurrent queries and keeps the
;;; ansync web service reasonably happy (though mysql is the bottleneck...
;;;
(def mysql-ds
     (dbcp/db-spec
      (let [ds (dbcp/mysql-datasource
                "127.0.0.1:3306" "refseq58" "root" @db-pw)]
        (dbcp/set-max-active! ds 50)
        (dbcp/set-min-max-idle! ds 5 20)
        ds)))


(defn create-ancestor-table []
  (sql/with-connection mysql-ds
    (sql/transaction
     (sql/create-table
      :ancestor
      [:ncbi_taxon_id :integer "PRIMARY KEY"]
      [:ancestors "varchar(1024)"]))))


(defn drop-ancestor-table []
  (sql/with-connection mysql-ds
    (sql/transaction
     (sql/drop-table :ancestor))))


(defn insert-ancestor-rows [rows]
  (sql/with-connection mysql-ds
    (sql/transaction
     (apply sql/insert-rows :ancestor rows))))



(defn populate-ancestor-table [file]
  (let [lazy-file (io/read-lines (io/file-str file))]
    (loop [v lazy-file
           cnt 0]
      (let [rows (map #(str/split #"\t\|\t" %1) (take 50 v))]
        (when (seq rows)
          ;;(prn rows)
          (insert-ancestor-rows rows)
          ;;(when (< cnt 500)
            (recur (drop 50 v)
                   (+ cnt 50)))))))



;;; (drop-ancestor-table)
;;; (create-ancestor-table)
;;; (populate-ancestor-table "/data2/BioData/taxdata/nodes-ancestor-names.txt")
