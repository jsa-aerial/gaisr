;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                        R D B . N E W - R N A S                           ;;
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

(ns edu.bc.bio.rdb.new-rnas
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
        [edu.bc.bio.sequtils.files
         :only [entry-parts read-seqs]]
        [edu.bc.log4clj
         :only [create-loggers log>]]
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        [edu.bc.bio.gaisr.db-actions
         :only [mysql-ds sql-query]])

  (:import javax.sql.DataSource
           com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource))




;;; ------------------------------------------------------------------------;;;
;;; Put the stuff into the DB ....

(defn get-gname-bioid-map []
  (let [snapshot "/data2/BioData/Archives/genome-id-map.clj"]
    (if (fs/exists? snapshot)
      (io/with-in-reader snapshot (read))
      (let [fields ["name" "bioentry_id"]
            where "name regexp \"^NC\""
            qmap (sql-query (str "select " (str/join ", " fields)
                                 " from bioentry "
                                 "where " where))
            gn-bid-map (into {} (map #(do [(% :name) (% :bioentry_id)])
                                     qmap))]
        (io/with-out-writer snapshot (prn gn-bid-map))
        gn-bid-map))))

(def *name-id-map* (get-gname-bioid-map))


(defmacro sql-update
  [sql-op]
  (let [sqlop? ((ns-interns 'clojure.contrib.sql) (first sql-op))
        opsym (when sqlop? (symbol "sql"  (name (first sql-op))))
        sql-op (if sqlop? (cons opsym (rest sql-op)) sql-op)]
    `(sql/with-connection mysql-ds
       (sql/transaction
        ~sql-op))))

(defn create-rnahit-tables []
  (sql-update
   (do-commands
    "CREATE TABLE rnahit (
          rnahit_id   INT(10) UNSIGNED NOT NULL auto_increment,
          name        varchar(40),
          gene_name   varchar(40),
          PRIMARY KEY (rnahit_id),
          UNIQUE (name)
     ) ENGINE=INNODB"
    "CREATE INDEX rnahit_gene ON rnahit(gene_name)"
    "CREATE INDEX rnahit_gene_name ON rnahit(name,gene_name)"

    "CREATE TABLE genome_rnahit (
        genome_rnahit_id   INT(10) UNSIGNED NOT NULL auto_increment,
        bioentry_id        INT(10) UNSIGNED NOT NULL,
        rnahit_id          INT(10) UNSIGNED NOT NULL,
        PRIMARY KEY (genome_rnahit_id)
     ) ENGINE=INNODB"
    "CREATE INDEX genome_rnahit_beid ON genome_rnahit(bioentry_id)"
    "CREATE INDEX genome_rnahit_rhid  ON genome_rnahit(rnahit_id)"
    )))

(defn drop-rnahit-tables []
  (sql-update
   (do-commands
    "drop table rnahit"
    "drop table genome_rnahit")))


(defn- insert-rnahits
  [values]
  (sql-update
   (apply
    sql/insert-values
    :rnahit
    [:name :gene_name] values)))

(defn update-rnahit
  [stospec]
  (let [stos (if (fs/directory? stospec)
               (sort (fs/directory-files stospec ".sto"))
               (list (fs/fullpath stospec)))]
    (insert-rnahits
     (partition-all
      2 (flatten
         (map #(let [rna-name (->> % fs/basename (str/split #"\.")
                                   first str/lower-case)
                     genes (->> % io/read-lines (drop 2) first
                                (str/split #"\s+") last
                                (str/split #" *, *"))]
                 (map (fn[g] [rna-name g])
                      genes))
              stos))))))

;;;(update-rnahit "/home/kaila/Bio/Tests/JSA/NewRNAs")
;;;(update-rnahit "/home/kaila/Bio/Tests/JSA/NewRNAs/RNA_00003.sto")

(defn- insert-genome-rnahits
  [values]
  (sql-update
   (apply
    sql/insert-values
    :genome_rnahit
    [:bioentry_id :rnahit_id] values)))

(defn update-genome-rnahit
  [stospec]
  (let [stos (if (fs/directory? stospec)
               (sort (fs/directory-files stospec ".sto"))
               (list (fs/fullpath stospec)))
        rnahit-ids (map #(->> % fs/basename
                              (str/split #"\.")
                              first str/lower-case
                              (pr-str)
                              (str "select rnahit_id from rnahit where name=")
                              (sql-query)
                              first (:rnahit_id))
                        stos)
        genome-groups (map
                       (fn[sto]
                         (->> sto (#(read-seqs % :info :name))
                              (map #(first (entry-parts %)))
                              (map #(*name-id-map* %))))
                       stos)]
    (map insert-genome-rnahits
         (map (fn[rnahit-id genomes]
                (map #(do [% rnahit-id]) genomes))
              rnahit-ids
              genome-groups))))

#_(update-genome-rnahit "/home/kaila/Bio/Tests/JSA/NewRNAs")

#_(map :name (sql-query "select name from rnahit"))

(def new-rna-genome-counts
     "select count(*) from
      (select distinct tx.taxon_id
              from bioentry as be,
                   taxon as tx,
                   ancestor as an,
                   rnahit as rh,
                   genome_rnahit as grh
              where be.bioentry_id=grh.bioentry_id
              and be.taxon_id=tx.taxon_id
              and be.name regexp \"^NC\"
              and be.description not regexp \"plasmid\"
              and tx.ncbi_taxon_id=an.ncbi_taxon_id
              and grh.rnahit_id=rh.rnahit_id
              and an.ancestors regexp \"/taxon/\"")

(def genome-counts
     "select count(*) from
      (select distinct tx.taxon_id
              from bioentry as be,
                   taxon as tx,
                   ancestor as an
              where be.taxon_id=tx.taxon_id
              and tx.ncbi_taxon_id=an.ncbi_taxon_id
              and be.name regexp \"^NC\"
              and be.description not regexp \"plasmid\"
              and an.ancestors regexp \"/taxon/\"")


(defn count-genomes
  [taxon & {:keys [new-rnas gene] :or {new-rnas false gene nil}}]
  (let [new-rna-pred " and rh.name="
        new-rna-clause (when (string? new-rnas)
                         (str new-rna-pred \" new-rnas \"))
        gene-clause (when gene
                      (str " and rh.gene_name=\"" gene \"))
        stmt (if new-rnas new-rna-genome-counts genome-counts)
        stmt (if new-rna-clause (str stmt new-rna-clause) stmt)
        stmt (if gene (str stmt gene-clause) stmt)
        stmt (str/replace-re #"/taxon/" taxon stmt)
        stmt (str stmt ") as foo")]
    (sql-query stmt)))


(comment
  (count-genomes "Firmicutes")
  (count-genomes "Bacteria")
  (count-genomes "Vibrionales")
  (count-genomes "Enterobacteriales")
  (count-genomes "endosymbiont")

  (count-genomes "Firmicutes" :new-rnas true)
  (count-genomes "Bacteria" :new-rnas true)
  (count-genomes "Vibrionales" :new-rnas true)
  (count-genomes "Enterobacteriales" :new-rnas true)
  (count-genomes "endosymbiont" :new-rnas true)
  )



(def bacterial-taxons-of-note
     ["Gammaproteobacteria" "Enterobacteriales" "Pseudomonadales"
      "Vibrionales" "Xanthomonadales" "Thiotrichales" "Oceanospirillales"
      "Pasteurellales" "Legionellales" "Alteromonadales" "Methylococcales"
      "Chromatiales" "Aeromonadales" "Legionellales" "Epsilonproteobacteria"
      "Deltaproteobacteria" "Alphaproteobacteria" "Betaproteobacteria"
      "Firmicutes" "Bacilli" "Clostridia" "Negativicutes" "Erysipelotrichi"
      "Spirochaetes" "Actinobacteria" "Verrucomicrobia" "Chlamydiae"
      "Cyanobacteria" "Bacteroidetes" "Chlorobi" "Fusobacteria"
      "Tenericutes" "Chloroflexi" "Thermotogae" "Deinococci" "Aquificae"
      "Acidobacteria" "Archaea" "Buchnera" "Blochmannia"
      "Wigglesworthia" "endosymbiont" "Xenorhabdus"])

(def bacterial-taxons-of-note-other
     ["Oceanospirillales" "Buchnera" "Blochmannia"
      "Wigglesworthia" "endosymbiont" "Xenorhabdus"])




#_(doseq [[gn cnt] (map #(do [% (-> % count-genomes  first vals first)])
                        bacterial-taxons-of-note)]
    (println (str gn ", " cnt)))


#_(doseq [[gn cnt] (map #(do [% (-> % (count-genomes :new-rnas true)
                                    first vals first)])
                        bacterial-taxons-of-note)]
    (println (str gn ", " cnt)))


#_(doseq [grp (for [rna (map :name (sql-query "select name from rnahit"))
                    taxon bacterial-taxons-of-note
                    :let [fullcnt (-> taxon count-genomes first vals first)
                          cnt (-> taxon (count-genomes :new-rnas rna)
                                  first vals first)]
                    :when (and (> cnt 0) (> fullcnt 0))]
                [rna taxon cnt (double (* 100.0 (/ cnt fullcnt)))])]
    (println (str/join ", " grp)))



#_(sort-by
   second
   (map
    #(do [(% :name) (last (str/split #"\s*,\s*" (% :ancestors)))])
    (sql-query
      "select be.name,an.ancestors
              from bioentry as be, taxon as tx, ancestor as an
              where be.taxon_id=tx.taxon_id
              and   tx.ncbi_taxon_id=an.ncbi_taxon_id
              and   an.ancestors regexp \"vibrionales\"
              and   be.description not regexp \"plasmid\"
              and   be.name regexp \"^NC\"")))

#_(sql-query
   "select * from bioentry where name=\"NC_000925\"")
