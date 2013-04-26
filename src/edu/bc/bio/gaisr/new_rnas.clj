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

(ns edu.bc.bio.gaisr.new-rnas
  (:require [clojure.contrib.sql :as sql]
            [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        [edu.bc.bio.sequtils.files
         :only [entry-parts read-seqs]]
        [edu.bc.log4clj
         :only [create-loggers log>]]
        [clojure.pprint
         :only [cl-format]]
        [edu.bc.bio.gaisr.db-actions
         :only [mysql-ds sql-query insts-by-rank]])

  (:import javax.sql.DataSource
           com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource))




;;; ------------------------------------------------------------------------;;;
;;; Put the stuff into the DB ....

(defn get-gname-bioid-map
  "Fetch a copy of genome names to bioentry ids.  Snapshot to file for
   future reference if snapshot does not exist.  Otherwise, read the
   snapshot at startup.
  "
  []
  (let [snapshot "/data2/BioData/Archives/genome-id-map.clj"]
    (if (fs/exists? snapshot)
      (io/with-in-reader snapshot (read))
      ;; Else query DB and cache results in snapshot
      (let [fields ["name" "bioentry_id"]
            where "name regexp \"^NC\""
            qmap (sql-query (str "select " (str/join ", " fields)
                                 " from bioentry "
                                 "where " where))
            gn-bid-map (into {} (map #(do [(% :name) (% :bioentry_id)])
                                     qmap))]
        (binding [*print-length* nil]
          (io/with-out-writer snapshot (prn gn-bid-map)))
        gn-bid-map))))

(defparameter "Map of genome names to bioentry ids"
  *name-id-map* (get-gname-bioid-map))




(defmacro sql-update
  [sql-op]
  (let [sqlop? ((ns-interns 'clojure.contrib.sql) (first sql-op))
        opsym (when sqlop? (symbol "sql"  (name (first sql-op))))
        sql-op (if sqlop? (cons opsym (rest sql-op)) sql-op)]
    `(sql/with-connection mysql-ds
       (sql/transaction
        ~sql-op))))


;;; CONSTRAINT unique_name UNIQUE (name,version)
(defn create-rnahit-table []
  (sql-update
   (do-commands
    "CREATE TABLE rnahit (
          rnahit_id     INT(10) UNSIGNED NOT NULL auto_increment,
          name          varchar(40) NOT NULL,
          version       INT(10) UNSIGNED NOT NULL,
          gene_name     varchar(40) NOT NULL,
          ncbi_taxon_id INT(11) UNSIGNED NOT NULL,
          PRIMARY KEY (rnahit_id)
     ) ENGINE=INNODB"
    "CREATE INDEX rnahit_gene ON rnahit(gene_name)"
    "CREATE INDEX rnahit_ncbi_tx_id ON rnahit(ncbi_taxon_id)"
    "CREATE INDEX rnahit_txs_genes ON rnahit(ncbi_taxon_id,gene_name)"
    "CREATE INDEX rnahit_gene_name ON rnahit(name,gene_name)"
    "CREATE INDEX rnahit_name_version ON rnahit(name,version)"
    )))

(defn create-verified-rnahit-table []
  (sql-update
   (do-commands
    "CREATE TABLE verified_rnahit (
          verified_id        INT(10) UNSIGNED NOT NULL auto_increment,
          rnahit_id          INT(10) UNSIGNED NOT NULL,
          verified_taxon_id  INT(11) UNSIGNED NOT NULL,
          verified_loc       varchar (20) NOT NULL,
          PRIMARY KEY (verified_id)
     ) ENGINE=INNODB"
    "CREATE INDEX verified_rhid     ON verified_rnahit(rnahit_id)"
    "CREATE INDEX verified_taxon_id ON verified_rnahit(verified_taxon_id)"
    "CREATE INDEX verified_location ON verified_rnahit(verified_loc)"
    "CREATE INDEX verified_tx_loc
            ON verified_rnahit(verified_taxon_id,verified_loc)"
    )))

(defn create-genome-rnahit-table []
  (sql-update
   (do-commands
    "CREATE TABLE genome_rnahit (
        genome_rnahit_id   INT(10) UNSIGNED NOT NULL auto_increment,
        bioentry_id        INT(10) UNSIGNED NOT NULL,
        rnahit_id          INT(10) UNSIGNED NOT NULL,
        PRIMARY KEY (genome_rnahit_id)
     ) ENGINE=INNODB"
    "CREATE INDEX genome_rnahit_beid ON genome_rnahit(bioentry_id)"
    "CREATE INDEX genome_rnahit_rhid  ON genome_rnahit(rnahit_id)"
    )))

(defn create-newhit-tables []
  (create-rnahit-table)
  (create-verified-rnahit-table)
  (create-genome-rnahit-table))


(defn drop-rnahit-table []
  (sql-update (do-commands "drop table rnahit")))

(defn drop-verified-rnahit-table []
  (sql-update (do-commands "drop table verified_rnahit")))

(defn drop-genome-rnahit-table []
  (sql-update (do-commands "drop table genome_rnahit")))

(defn drop-newhit-tables []
  (drop-rnahit-table)
  (drop-verified-rnahit-table)
  (drop-genome-rnahit-table))


;;; (do
;;;   (drop-rnahit-table)
;;;   (create-rnahit-table)
;;;   (drop-verified-rnahit-table)
;;;   (create-verified-rnahit-table))


(defn get-sto-files [designator]
  (cond
   (coll? designator) (map fs/fullpath designator)
   (fs/directory? designator) (sort (fs/directory-files designator ".sto"))
   :else (list (fs/fullpath designator))))

(defn get-mlab-gfinfo [sto]
  (let [rna-name (->> sto fs/basename (str/split #"\.")
                      first str/lower-case)
        [rna-name version] (str/split #"-" rna-name)
        gfs (take-until
             #(not= "#=GF " (str/take 5 %))
             (drop 1 (io/read-lines sto)))
        infomap  (into
                  {}
                  (for [gf gfs
                        :let [bits (str/split #"\s+" gf)
                              kind (->> bits second str/lower-case keyword)
                              data (third bits)]
                        :when (in kind [:rna :taxon :in-vivo :in-vitro])]
                    [kind (cond
                           (= kind :rna) data
                           (= kind :taxon) (Integer. data)
                           :else (map #(Integer. %)
                                      (when data (str/split #", *" data))))]))
        rna (infomap :rna)
        taxon (infomap :taxon)
        invitros (infomap :in-vitro)
        invivos (infomap :in-vivo)]
    (if (and (empty? invitros) (empty invivos))
      [[rna-name version rna taxon 0 "na"]]
      (reduce (fn[v vid]
                (conj v [rna-name version rna taxon vid "vitro"]))
              (reduce (fn[v vid]
                        (conj v [rna-name version rna taxon vid "vivo"]))
                      [] invivos)
              invitros))))

;;;(map get-mlab-gfinfo (get-sto-files "/home/kaila/Bio/Tests/JSA/NewRNAs"))
;;;(map get-mlab-gfinfo (get-sto-files "/home/kaila/Bio/Tests/JSA"))


(defn- insert-rnahits
  [values]
  (sql-update
   (apply
    sql/insert-values
    :rnahit
    [:name :version :gene_name :ncbi_taxon_id] values)))

(defn update-rnahit
  [stospec]
  (let [stos (get-sto-files stospec)]
    (insert-rnahits
     (sort-by first (set (map (fn[[n v gn ntx & tail]] [n v gn ntx])
                              (apply concat (map get-mlab-gfinfo stos))))))))

;;;(update-rnahit "/home/kaila/Bio/Tests/JSA/NewRNAs")
;;;(update-rnahit "/home/kaila/Bio/Tests/JSA")
;;;(update-rnahit "/home/kaila/Bio/Tests/JSA/V2")
;;;(update-rnahit "/home/kaila/Bio/Tests/JSA/NewRNAs/RNA_00003.sto")
;;;
;;;(update-rnahit "/home/kaila/Bio/Tests/RNA_00011-3.sto")

(def rnaname-id-map
     (reduce (fn[M m]
               (assoc M (str (m :name) "-" (m :version))
                      (m :rnahit_id)))
             {} (sql-query
                 "select rh.name,rh.version,rh.rnahit_id from rnahit as rh")))

(defn- insert-verified-info
  [values]
  (sql-update
   (apply
    sql/insert-values
    :verified_rnahit
    [:rnahit_id :verified_taxon_id :verified_loc] values)))

(defn update-verified-info
  [stospec]
  (let [stos (get-sto-files stospec)]
    (insert-verified-info
     (sort-by first (map (fn[[n v gn ntx vid vloc]]
                           [(rnaname-id-map (str n "-" v)) vid vloc])
                         (apply concat (map get-mlab-gfinfo stos)))))))

;;;(update-verified-info "/home/kaila/Bio/Tests/JSA/NewRNAs")
;;;(update-verified-info "/home/kaila/Bio/Tests/JSA")
;;;(update-verified-info "/home/kaila/Bio/Tests/JSA/V2/")
;;;(update-verified-info "/home/kaila/Bio/Tests/RNA_00011-3.sto")


(defn- insert-genome-rnahits
  [values]
  (sql-update
   (apply
    sql/insert-values
    :genome_rnahit
    [:bioentry_id :rnahit_id] values)))

(defn update-genome-rnahit
  [stospec]
  (let [stos (get-sto-files stospec)
        rnahit-ids (map #(->> % fs/basename
                              (str/split #"\.")
                              first str/lower-case
                              rnaname-id-map)
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

;;; (update-genome-rnahit "/home/kaila/Bio/Tests/JSA/NewRNAs")
;;; (update-genome-rnahit "/home/kaila/Bio/Tests/JSA/NewRNAs/RNA_00006.sto")
;;; (update-genome-rnahit "/home/kaila/Bio/Tests/JSA")
;;; (update-genome-rnahit "/home/kaila/Bio/Tests/JSA/V2/RNA_00012-2.sto")
;;; (update-genome-rnahit "/home/kaila/Bio/Tests/RNA_00011-3.sto")


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
  [taxon & {:keys [new-rnas version gene]
            :or {new-rnas nil version nil gene nil}}]
  (let [new-rna-name (when new-rnas (str " and rh.name=\"" new-rnas \"))
        new-rna-version (when version (str " and rh.version=" version))
        new-rna-clause (when (or new-rna-name new-rna-version)
                         (str new-rna-name new-rna-version))
        gene-clause (when gene
                      (str " and rh.gene_name=\"" gene \"))
        stmt (if new-rnas new-rna-genome-counts genome-counts)
        stmt (if new-rna-clause (str stmt new-rna-clause) stmt)
        stmt (if gene (str stmt gene-clause) stmt)
        stmt (str/replace-re #"/taxon/" taxon stmt)
        stmt (str stmt ") as foo")]
    ;;(println stmt)
    (sql-query stmt)))


(def new-rna-taxon-genome-subtbl
     "(select distinct tx.taxon_id
              from bioentry as be,
                   taxon as tx,
                   ancestor as an,
                   rnahit as rh,
                   genome_rnahit as grh
              where be.bioentry_id=grh.bioentry_id
              and be.taxon_id=tx.taxon_id
              and tx.ncbi_taxon_id=an.ncbi_taxon_id
              and grh.rnahit_id=rh.rnahit_id
              and an.ancestors regexp \"/taxon/\"")

(def taxon-genome-subtbl
     "(select distinct tx.taxon_id
              from bioentry as be,
                   taxon as tx,
                   ancestor as an
              where be.taxon_id=tx.taxon_id
              and tx.ncbi_taxon_id=an.ncbi_taxon_id
              and an.ancestors regexp \"/taxon/\"")


(defn get-genomes
  [taxon & {:keys [plasmids ncs-only new-rnas version gene]
            :or {plasmids false ncs-only true
                 new-rnas nil version nil gene nil}}]
  (let [plasmid-clause " and be.description not regexp \"plasmid\""
        ncsonly-clause " and be.name regexp \"^NC\""
        rna-name       (when new-rnas (str " and rh.name=\"" new-rnas \"))
        rna-version    (when version (str " and rh.version=" version))
        rna-clause     (when (or rna-name rna-version)
                         (str rna-name rna-version))
        gene-clause    " and rh.gene_name="
        gene-clause (when gene (str gene-clause  \" gene \"))

        subtbl (if rna-clause new-rna-taxon-genome-subtbl taxon-genome-subtbl)
        subtbl (str/replace-re #"/taxon/" taxon taxon-genome-subtbl)
        subtbl (if plasmids subtbl (str subtbl plasmid-clause))
        subtbl (if ncs-only (str subtbl ncsonly-clause) subtbl)
        subtbl (str subtbl ") as foo ")

        tables " from bioentry as be,"
        tables (str (if (or rna-clause gene-clause)
                      (str tables "rnahit as rh,genome_rnahit as grh,")
                      tables)
                    subtbl)
        preds " where be.taxon_id=foo.taxon_id"
        preds (if (or rna-clause gene-clause)
                (str preds
                     " and be.bioentry_id=grh.bioentry_id"
                     " and grh.rnahit_id=rh.rnahit_id")
                preds)

        stmt (str "select be.name, be.bioentry_id, foo.taxon_id" tables preds)
        stmt (if ncs-only (str stmt ncsonly-clause) stmt)
        stmt (if plasmids stmt (str stmt plasmid-clause))
        stmt (if rna-clause (str stmt rna-clause) stmt)
        stmt (if gene-clause (str stmt gene-clause) stmt)]
    (sql-query stmt)))


(defn get-ncs-by-taxon-rna-gene
  "For TAXON (a single string naming a taxon or a collection of such
   strings) and for the options RNA and/or GENE return a collection
   comprised of [tx rn ncs] or [tx gn ncs] or [tx rn gn ncs], where:

   * Collection contains only one of these element types

   * tx is one of TAXON

   * rn is a new rna name (eg, \"rna_00001\")

   * gn is a gene name (eg, \"rplA\"

   * ncs is a collection of NC* genome names
  "

  [taxon & {:keys [rna gene] :or {rna nil gene nil}}]
  {:pre [(or (string? rna) (string? gene))]}
  (for [taxon (ensure-vec taxon)
        :let [stmt "select be.name,rh.name as rna, rh.gene_name as gene
                           from bioentry as be,
                                ancestor as an,
                                taxon as tx,
                                rnahit as rh,
                                genome_rnahit as grh "
              stmt (str stmt "where "
                             "rh.rnahit_id=grh.rnahit_id "
                             (if rna (str "and rh.name=\"" rna "\" ") "")
                             (if gene (str "and rh.gene_name=\"" gene "\" ") "")
                             "and   grh.bioentry_id=be.bioentry_id
                              and   be.taxon_id=tx.taxon_id
                              and   tx.ncbi_taxon_id=an.ncbi_taxon_id
                              and   an.ancestors regexp \"")
              stmt (str stmt taxon "\"")
              ;;_ (println stmt)
              result (sql-query stmt)
              ncs (map :name result)
              rnas (map :rna result)
              rnas-ncs (reduce (fn[m e]
                                 (let [k (:rna e)
                                       v (:name e)]
                                   (assoc m k (conj (get m k []) v))))
                               {} result)
              genes-ncs (reduce (fn[m e]
                                 (let [k (:gene e)
                                       v (:name e)]
                                   (assoc m k (conj (get m k []) v))))
                               {} result)]]
    (cond
     (and rna gene) [taxon rna gene ncs]
     gene [taxon gene rnas-ncs]
     rna [taxon rna genes-ncs])))


(defn print-ncs-by-taxon
  [coll]
  (if (= (count (first coll)) 4)
    (doseq [[taxon rn gn ncs] coll]
      (println (str/join ", " [taxon rn gn]) ":\n")
      (doseq [nc ncs] (println (str "    " nc))))
    (doseq [[taxon x ncs] coll]
      (println "\n" (str/join ", " [taxon x]) ":\n")
      (doseq [[n vs] ncs]
        (println (str "    " n " -> "))
        (doseq [v vs] (println (str "        " v)))))))


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

(def sorted-bacterial-taxons-of-note
     ["Acidobacteria" "Actinobacteria" "Aeromonadales" "Alphaproteobacteria"
      "Alteromonadales" "Aquificae" "Archaea" "Bacilli" "Bacteroidetes"
      "Betaproteobacteria" "Blochmannia" "Buchnera" "Chlamydiae" "Chlorobi"
      "Chloroflexi" "Chromatiales" "Clostridia" "Cyanobacteria" "Deinococci"
      "Deltaproteobacteria" "Enterobacteriales" "Epsilonproteobacteria"
      "Erysipelotrichi" "Firmicutes" "Fusobacteria" "Gammaproteobacteria"
      "Legionellales" "Legionellales" "Methylococcales" "Negativicutes"
      "Oceanospirillales" "Pasteurellales" "Pseudomonadales" "Spirochaetes"
      "Tenericutes" "Thermolithobacteria" "Thermotogae" "Thiotrichales"
      "Verrucomicrobia" "Vibrionales" "Wigglesworthia" "Xanthomonadales"
      "Xenorhabdus" "endosymbiont"])




;;;  rnas ["rna_00011" "rna_00012" "rna_00013"]] ;(keys rnaname-id-map)
;;;
(defn rna-taxon-info
  "For each rna in RNAS, a single, or collection of, new rna
   name-version items (eg RNA_00011-2) create a breakdown of count,
   percentage of total and NC list for each taxon in TAXONS.  Place
   results in outfile.
  "
  [out-file & {:keys [rnas taxons]
               :or {rnas (sort (keys rnaname-id-map))
                    taxons sorted-bacterial-taxons-of-note}}]
  (io/with-out-writer out-file
    (doseq [grp (for [rna-v (ensure-vec rnas)
                      taxon taxons
                      :let [[rna v] (str/split #"-" rna-v)
                            fullcnt (-> taxon count-genomes first vals first)
                            cnt (-> taxon
                                    (count-genomes :new-rnas rna :version v)
                                    first vals first)]
                      :when (and (> cnt 0) (> fullcnt 0))]
                  [rna (str "v" v) taxon cnt
                   fullcnt (double (* 100.0 (/ cnt fullcnt)))
                   (->> (get-ncs-by-taxon-rna-gene taxon :rna rna)
                        first third vals first set)])]
      (println)
      (println (str/join ", " (take 6 grp)))
      (doseq [nc (last grp)]
        (println nc))))
  out-file)


(comment
  (rna-taxon-info "/home/kaila/Bio/Tests/nc-12-taxinfo.txt"
                  :rnas "rna_00012-2")
  (str/split #"-" "rna_00012-2")
  (count-genomes "Acidobacteria" :new-rnas "rna_00012" :version 2)

  (rna-taxon-info "/home/kaila/Bio/Tests/nc-12-2-and-13-1-taxinfo.txt"
                  :rnas ["rna_00012-2" "rna_00013-1"]
                  :taxons (conj sorted-bacterial-taxons-of-note
                                "Thermoanaerobacterales"
                                "Natranaerobiales"
                                "Clostridiales"
                                "Halanaerobiales"
                                "Bacillales"
                                "Lactobacillales"))
  )


(comment
  (sto-intersect "/home/kaila/Bio/Tests/JSA/NewRNAs/RNA_00006.sto"
                 "/home/kaila/Bio/Tests/JSA/RNA_00012.sto")

  (sto-intersect "/home/kaila/Bio/Tests/JSA/NewRNAs/RNA_00010.sto"
                 "/home/kaila/Bio/Tests/JSA/RNA_00013.sto")

  (count-genomes "Firmicutes")
  (count-genomes "Bacteria")
  (count-genomes "Vibrionales")
  (count-genomes "Enterobacteriales")
  (count-genomes "endosymbiont")
  (get-genomes "endosymbiont")
  (count (get-genomes "Firmicutes"))
  (count (get-genomes "Bacteria"))

  (count (get-genomes "Firmicutes" :new-rnas "rna_00013"))
  (count (->> (get-ncs-by-taxon-rna-gene "Firmicutes" :rna "rna_00013")
              first third (#(% "rplT"))))
  (count-genomes "Bacteria" :new-rnas true)
  (count-genomes "Vibrionales" :new-rnas true)
  (count-genomes "Enterobacteriales" :new-rnas true)
  (count-genomes "endosymbiont" :new-rnas true)
  )




#_(doseq [[gn cnt] (map #(do [% (-> % count-genomes  first vals first)])
                        sorted-bacterial-taxons-of-note)]
    (println (str gn ", " cnt)))


#_(doseq [[gn cnt] (map #(do [% (-> % (count-genomes :new-rnas true)
                                    first vals first)])
                        sorted-bacterial-taxons-of-note)]
    (println (str gn ", " cnt)))



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


#_(sql-query
   "select be.name,an.ancestors
           from bioentry as be, taxon as tx, ancestor as an
           where be.taxon_id=tx.taxon_id
           and   tx.ncbi_taxon_id=an.ncbi_taxon_id
           and   an.ancestors regexp \"Leuconostoc\"
           and   be.name regexp \"^NC\"")

#_(sql-query
   "select be.name,an.ancestors
           from bioentry as be, taxon as tx, ancestor as an
           where be.taxon_id=tx.taxon_id
           and   tx.ncbi_taxon_id=an.ncbi_taxon_id
           and   an.ancestors regexp \"Firmicute\"
           and   be.name regexp \"^NC\"")


#_(print-ncs-by-taxon
   (get-ncs-by-taxon-rna-gene
    ["Deltaproteobacteria" "Epsilonproteobacteria"
     "Alphaproteobacteria" "Betaproteobacteria"]
    :rna "rna_00007"))

#_(print-ncs-by-taxon
   (get-ncs-by-taxon-rna-gene
    ["Tenericutes" "Thermotogae" "Fusobacteria" "Chloroflexi"
     "Cyanobacteria" "Verrucomicrobia/Chlamydia" "Deinococci"
     "Bacteroidetes" "Chlorobi"]
    :gene "rplA"))

#_(let [ncs (-> (get-ncs-by-taxon-rna-gene "Firmicutes" :gene "rplA")
                first third
                (get "rna_00007")
                (#(reduce (fn[m nm] (assoc m nm nm)) {} %)))]
    (->> (read-seqs "/home/kaila/Bio/Tests/RNA_00007.sto" :info name)
         (map entry-parts)
         (map first)
         (filter #(get ncs %))
         (reduce (fn [m nm]
                   (assoc m nm (inc (get m nm 0))))
                 {})
         (group-by #(val %))
         (#(do [(map first (% 1)) (map first (% 2))]))
         ((fn[[ones twos]]
            (let [one-cnt (count ones)
                  two-cnt (count twos)
                  total (+ one-cnt two-cnt)]
              (/ two-cnt total))))))


#_(print-ncs-by-taxon
   (get-ncs-by-taxon-rna-gene
    " Thermus "
    :rna "rna_00009"))

#_(print-ncs-by-taxon
   (get-ncs-by-taxon-rna-gene
    "Tenericutes"
    :gene "rplA"))


(comment
  (insts-by-rank "Bacillus subtilis" "family")
  (count (insts-by-rank
          "Proteobacteria" "genus" :pred #(re-find #"^NC" (:name %))))
  )