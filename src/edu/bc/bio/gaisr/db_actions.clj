;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                            D B - A C T I O N S                           ;;
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

(ns edu.bc.bio.gaisr.db-actions
  (:require [clojure.contrib.sql :as sql]
            [org.bituf.clj-dbcp :as dbcp]
            [clojure.contrib.string :as str]
            [clojure.set :as set]
            [clojure.contrib.json :as json]
            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        [edu.bc.log4clj :only [create-loggers log>]]

        [clojure.pprint
         :only [cl-format]])

  (:import javax.sql.DataSource
           com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource
           (java.nio ByteBuffer CharBuffer)
           (java.io PushbackReader InputStream InputStreamReader
                    FileInputStream)
           ))


(def db-pw (atom "rna314rulz"))

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


(defn sql-query [stmt & {:keys [f p] :or {f identity p false}}]
  (when p (println stmt))
  (sql/with-connection mysql-ds
    (sql/with-query-results qresults [stmt]
      (let [results (doall qresults)]
        (f results)))))




;;; -------------------- Operon Location Query ----------------------;;;


(defparameter *bogus-counts* (atom {}))

;;; For a given bioentry (organism), fetch the operon location
;;; information that falls between start and end.
;;;
(defn operon-location-query
  "For a given bioentry (organism) given by ENTRY-NAME (an NCBI genome name),
   fetch the operon location information that falls between START and END.
   START and END are optional.  START defaults to 1.  END defaults to large
   number exceeding any bacterial genome size (should default to
   Long/MAX_VALUE, but MySql chokes on this...)"
  ([entry-name]
     (operon-location-query entry-name 1 100000000))
  ([entry-name start]
     (operon-location-query entry-name start 100000000))
  ([entry-name start end]
     (let [start (str start)
           end (str end)
           fields (str "op.operon_id, "
                       "oploc.start_pos as start, "
                       "oploc.end_pos as end, "
                       "oploc.strand")
           tables (str "bioentry as be, "
                       "operon as op, "
                       "operon_loc as oploc")
           constraints [(str "be.name=" \" entry-name \")
                        "be.bioentry_id=op.bioentry_id"
                        "op.operon_id=oploc.operon_id"
                        (str "oploc.start_pos between " start " and " end)
                        (str "oploc.end_pos between " start " and " end)]
           options     "order by start"
           constraint  (str/join " and " constraints)
           stmt (str "select " fields
                     " from " tables
                     " where " constraint
                     " " options)]
       (sql/with-connection mysql-ds
         (sql/with-query-results qresults [stmt]
           (doall
            (filter
             #(if (> (- (% :end) (% :start)) 0)
                true
                (do (swap! *bogus-counts*
                           (fn[m]
                             (assoc m [entry-name (% :start) (% :end)] 1)))
                    nil))
             qresults)))))))




;;; -------------------- Genome context Feature Query ----------------------;;;


(def multival-key
     #{:strand :length :name})

(defn sift-features [k curval nval]
  (let [vals (ensure-vec curval)]
    (cond
     (multival-key k) (conj vals nval)
     (in nval vals) vals
     true (conj vals nval))))

(defn merge-feature-kvs [feature-group]
  (apply merge-with* sift-features (second feature-group)))

(defn combine-vals2name [names vals]
  (reduce (fn[m [n v]]
            (assoc m n (let [x (get m n)] (if x (str x ", " v) v))))
          {}
          (map #(do [%1 %2]) names vals)))

(defn make-feature-struct [fmap]
  (let [sfid    (first (ensure-vec (fmap :seqfeature_id)))
        sftype  (first (ensure-vec (fmap :sftype)))
        starts  (ensure-vec (fmap :start_pos))
        ends    (ensure-vec (fmap :end_pos))
        strands (ensure-vec (fmap :strand))
        locs    (map (fn [s e r]
                       {:start s :end e
                        :len (inc (- e s))
                        :strand r})
                     starts ends strands)
        names   (ensure-vec (fmap :name))
        vals    (ensure-vec (fmap :value))
        nvmap   (combine-vals2name names vals)
        nvs     (map #(hash-map :name (key %1) :value (val %1)) nvmap)
        nvs     (sort-by #(%1 :name) nvs)]
    {:sfid sfid :sftype sftype :locs locs :nvs nvs}))



(defn feature-query [args]
  (let [name (args :name)
        id   (args :id)
        gbid (args :gbid)]
    (log> "rootLogger" :info
          "Feature Query: ~A:~A" name id)
    (let [fields (str "sf.seqfeature_id,sfterm.name as sftype, "
                      "start_pos, end_pos, end_pos-start_pos+1 as length, "
                      "strand, term.name, value")
          tables (str "seqfeature as sf, "
                      "(select term_id, name from term) as sfterm, "
                      "location as loc, "
                      "seqfeature_qualifier_value as sfqv, "
                      "term")
          constraints [(str "sf.bioentry_id=" id)
                       "sf.type_term_id=sfterm.term_id"
                       "sf.seqfeature_id=loc.seqfeature_id"
                       "sf.seqfeature_id=sfqv.seqfeature_id"
                       "sfqv.term_id=term.term_id"
                       "term.term_id not in (255, 256, 262)"]
          options     "order by sf.seqfeature_id"
          constraint (str/join " and " constraints)
          stmt (str "select " fields
                    " from " tables
                    " where " constraint
                    " " options)]
      (sql/with-connection mysql-ds
        (sql/with-query-results qresults [stmt]
          {:type :json
           :body (json/json-str
                  (conj ; First element is the gid, for ncbi queries
                   (->> (doall qresults)
                       (group-by #(% :seqfeature_id))
                       (sort (fn [l r] (< (first l) (first r))))
                       (map merge-feature-kvs)
                       (map make-feature-struct)
                       (sort-by #(%1 :sftype)))
                   gbid))}
          )))))




;;; --------------------- Hit Context Feature Query ----------------------- ;;;


(def +start-delta+ 5000)
(def +end-delta+ 5000)

(def seqfeature-subtable
     "(select sf.seqfeature_id,sf.type_term_id,term.name,
              loc.start_pos,loc.end_pos,loc.strand
           from bioentry as be,
                seqfeature as sf,
                location as loc,
                term
           where be.name=\"/xxx/\"
           and sf.bioentry_id=be.bioentry_id
           and sf.type_term_id=term.term_id
           and sf.seqfeature_id=loc.seqfeature_id
           and loc.start_pos between /spre/ and /epost/
           and loc.end_pos between /spre/ and /epost/)")


;;; For a given bioentry (genome) + "hit" start..end subseq location
;;; on its genome (as obtained by some "comparison" method to a
;;; candidate sequence), obtain the feature sets by location
;;; constraints.  The location constraints are (- start +start-delta+)
;;; and (+ end +end-delta+).  Any feature found within the deltas is
;;; obtained.
;;;
(defn hit-features-query [entry-name start end]
  (let [x (- start +start-delta+)
        spre (str (if (pos? x) x 1))
        epost (str (+ end +end-delta+))
        start (str start)
        end (str end)
        subtable (str/replace-re #"/xxx/" entry-name seqfeature-subtable)
        subtable (str/replace-re #"/spre/" spre subtable)
        subtable (str/replace-re #"/epost/" epost subtable)
        fields (str "sfqv.seqfeature_id, "
                    "x.type_term_id, x.name as sftype, "
                    "x.start_pos, x.end_pos, "
                    "x.end_pos-x.start_pos+1 as length, x.strand, "
                    "term.term_id, term.name, sfqv.value")
        tables (str "seqfeature_qualifier_value as sfqv, "
                    "term, "
                    subtable " as x")
        constraints ["x.seqfeature_id=sfqv.seqfeature_id"
                     "term.term_id=sfqv.term_id"]
        options     "order by sfqv.seqfeature_id"
        constraint  (str/join " and " constraints)
        stmt (str "select " fields
                  " from " tables
                  " where " constraint
                  " " options)]
    (sql/with-connection mysql-ds
      (sql/with-query-results qresults [stmt]
        (->> (doall qresults)
             (group-by #(% :seqfeature_id))
             (sort (fn [l r] (< (first l) (first r))))
             (map merge-feature-kvs)
             (map make-feature-struct)
             (sort-by #(%1 :sftype)))))))


(defn base-info-query [names]
  (let [fields (str  "be.name, be.version, be.description, "
                     "be.identifier as gbid, be.bioentry_id, "
                     "a.ncbi_taxon_Id as taxon_id, "
                     "tn.name as taxname, a.ancestors")
        tables "bioentry as be, taxon, taxon_name as tn, ancestor as a "
        constraints [(str "be.name in (" (str/join ", " names) ")")
                     "be.taxon_id=taxon.taxon_id"
                     "taxon.taxon_id=tn.taxon_id"
                     "tn.name_class=\"scientific name\""
                     "taxon.ncbi_taxon_id=a.ncbi_taxon_id"]
        constraint (str/join " and " constraints)
        stmt (str "select " fields
                  " from " tables
                  " where " constraint)]
    (sql/with-connection mysql-ds
      (sql/with-query-results qresults [stmt]
        (doall qresults)))))




;;; ----------------------- Taxonomic Based Querys ------------------------ ;;;


;;; This is a set holding the DB fields to return out of a query.
;;; Clojure vars, being Lisp dynamic vars, are dynamically rebindable
;;; per any call stack (including per thread and differnt points in
;;; thread call trees).  So, this is easily dynamically changed and
;;; changable per query type by rebinding in the handlers for the
;;; different queries.
;;;
(defparameter *ret-keys*
     (hash-set :taxname :taxon_id
               :bioentry_id :identifier :name :description
               :version :sfcount :ancestors))

;;; Take the result map from the taxon driven bioentries query, and
;;; transform the contents.  This revolves around keeping only the
;;; keys listed in the *RET-KEYS* map and further, also transforms the
;;; taxonomy ancestors result to filter out the top two nodes
;;; (cellular org and bacteria as they are _always_ part of _our_
;;; results) and the final node which is the originating taxon which
;;; is included as a separate field.
;;;
(defn xform-result [result-map]
  (into {} (keep
            (fn [kv]
              (let [v (*ret-keys* (key kv))]
                (when v
                  (cond
                   (= :ancestors v)
                   (let [ans (vec (str/split #", " (val kv)))
                         ans (subvec ans 2 (dec (count ans)))]
                     [(key kv) (str/join ", " ans)])

                   (= :identifier v) [:gbid (val kv)]

                   :else kv))))
            result-map)))


;;; Obtains the total count of "features" for bioentries.  The feature
;;; count for a bioentry is simply the count of corresponding rows
;;; with the bioentries id in the seqfeature table. This is BIOSQL
;;; SPECIFIC
;;;
(def seqfeature-count-table
     "(select bioentry_id, count(*) as sfcount
       from seqfeature group by bioentry_id)")


;;; Try to get a root taxon from taxonomy based on a name or partial
;;; name.  Selects minimum taxon-id from nodes whose name regexp
;;; matches the substitution term for xxx.
;;;
;;; Known Bugs: This can misfire if the nodes of the subtrees are not
;;; ordered by tree level by taxon_id.
;;;
(def root-taxon-table
     "(select min(taxon.taxon_id) as taxid
             from taxon_name as tn,
                  taxon
             where tn.name regexp \"^/xxx/\"
             and tn.name_class=\"scientific name\"
             and tn.taxon_id=taxon.taxon_id)")


;;; Sub query (table) for overall taxon scoped bioentries query.  This
;;; query attempts (and mostly succeeds) to obtain the set of all
;;; taxon categories under a given taxon category with something like
;;; "reasonable" performance.
;;;
;;; Since we are currently using an RDB (biosql) this is never going
;;; to be great.  Unless you did a full "closure" precompute of all
;;; paths - but that would result in a table with 100s of millions of
;;; rows and kill MySql outright.
;;;
;;; The search centers around nodes whose taxon name match the REGEX
;;; that is substituted for the XXX pattern and comes from the query
;;; input on the main web page.
;;;
;;; Known Bugs: This can misfire if the nodes associated with finding
;;; the root of the subtree we are after are not ordered by tree level
;;; by taxon_id.
;;;
(def taxon-nodes-table
    "(select tn.name, nodes.*
      from taxon_name as tn,
           taxon as parent,
           taxon as nodes,
           (select min(taxon.taxon_id) as taxid, name
                   from taxon_name as tn,
                        taxon
                   where tn.name regexp \"^/xxx/\"
                   and tn.name_class=\"scientific name\"
                   and tn.taxon_id=taxon.taxon_id) as root
      where tn.taxon_id=nodes.taxon_id
      and tn.name_class=\"scientific name\"
      and nodes.left_value between parent.left_value and parent.right_value
      and parent.taxon_id=root.taxid)")

(def taxon-nodes-query
     "select tn.name, nodes.*
      from taxon_name as tn,
           taxon as parent,
           taxon as nodes
      where tn.taxon_id=nodes.taxon_id
      and tn.name_class=\"scientific name\"
      and nodes.left_value between parent.left_value and parent.right_value
      and parent.taxon_id=xxx")

(def taxon-leafs-query
     "select taxon_name.name, leafs.*
       from taxon_name,
            taxon as parent,
            taxon as leafs
       where taxon_name.taxon_id=leafs.taxon_id
       and taxon_name.name_class=\"scientific name\"
       and leafs.left_value between parent.left_value and parent.right_value
       and leafs.left_value + 1 = leafs.right_value
       and parent.taxon_id=xxx")




(defn names->tax [names]
  (let [names (map #(cl-format nil "~S" %) (keys names))
        fields "be.name, a.ancestors"
        tables ["bioentry as be"
                "taxon as tx"
                "ancestor as a"]
        tables (str/join "," tables)
        constraints [(str "be.name in ("
                          (str/join "," names) ")")
                     "be.taxon_id=tx.taxon_id"
                     "tx.ncbi_taxon_id=a.ncbi_taxon_id"]
        constraint (str/join " and " constraints)
        stmt (str "select " fields
                  " from " tables
                  " where " constraint)]
    ;;(prn stmt)
    (sql/with-connection mysql-ds
      (sql/with-query-results qresults [stmt]
        (doall qresults)))))




(defn taxnode [term]
  (let [taxid (cond
               (integer? term) (str term)
               (intstg? term) term
               :else
               (sql-query
                (let [c1 (str/take 1 term)
                      term-pat (cl-format nil "[~A~A]~A"
                                          c1 (str/swap-case c1)
                                          (str/drop 1 term))]
                  (str/replace-re #"/xxx/" term-pat root-taxon-table))
                :f #(:taxid (first %))))
        fields "tn.name, nodes.*"
        tables "taxon as nodes, taxon_name as tn"
        constraints [(str "nodes.taxon_id=" taxid)
                     "tn.taxon_id=nodes.taxon_id"
                     "tn.name_class=\"scientific name\""]
        stmt (str "select " fields " from " tables
                  " where " (str/join " and " constraints))]
    (sql-query stmt :f first)))


(defn xform-taxnode [node]
  ;; sql spins off new thread so must bind *ret-keys* in that
  ;; thread and further, must force seq eval in that thread!
  (binding [*ret-keys* (conj *ret-keys* :ncbi_taxon_id :node_rank)]
    (doall (map xform-result node))))


(defn taxon-subclasses [term]
  (let [taxid (str (:taxon_id (taxnode term)))]
      (sql-query (str/replace-re #"xxx" taxid taxon-nodes-query)
                 :f xform-taxnode)))

(defn taxon-leafs [term]
  (let [taxid (str (:taxon_id (taxnode term)))]
    (sql-query (str/replace-re #"xxx" taxid taxon-leafs-query)
               :f xform-taxnode)))

(defn taxon-direct-instances [term]
  (let [taxid (str (:taxon_id (taxnode term)))
        stmt (str "select * from bioentry where taxon_id=" taxid)]
    (sql-query stmt)))


(defn insts-by-rank
  "Fetch the set of instances (genomes) classified by taxon by rank.

   TERM is a taxon denotation - an integer or integer string that is
   the taxon_id, or a name or partial name that is a
   prefix (potentially full) of the taxon's name.

   RANK is the name of a taxonomic rank under the taxon denoted by
   TERM.

   PRED is a predicate to filter the instance sets of each rank
   before random selection of instances from them.  Defaults to
   identity - no filtering.

   ENTRY-TYPE is one of #{:NC :NZ :OTHER :ALL} indicating the type of
   genome entry to return.  :NC indicates NC* genomes (annotated), :NZ
   indicating metagenomic (of sub type???), :OTHER is any non NC* or
   NZ*, and :ALL is any type.

   Result is a set of 'node' entries.  Each such entry is a map of
   information on the organism, including its bioentry id, name, taxon
   id, ncbi taxon id, et.al.
   "
  [term rank & {:keys [pred entry-type] :or {pred identity entry-type :NC}}]
  (keep (fn[s] (filter pred (seq s)))
        (map #(flatten
               (keep (fn[m]
                       (let [insts (taxon-direct-instances (m :taxon_id))]
                         (when insts insts)))
                     %))
             (map #(taxon-leafs (% :taxon_id))
                  (filter #(= (:node_rank %) rank)
                          (taxon-subclasses term))))))


(defn rand-insts-by-rank
  "Fetch a 'random' set of instances (genomes) classified by taxon by rank.

   TERM is a taxon denotation - an integer or integer string that is
   the taxon_id, or a name or partial name that is a
   prefix (potentially full) of the taxon's name.

   RANK is the name of a taxonomic rank under the taxon denoted by
   TERM.

   PRED is a predicate to filter the instance sets of each rank
   before random selection of instances from them.  Defaults to
   identity - no filtering.

   CNT is the size of the randomly selected subsets from each rank.

   Result is a set of 'node' entries formed by the union of the
   randomly selected subsets of each rank.  Each such entry is a map
   of information on the organism, including its bioentry id, name,
   taxon id, ncbi taxon id, et.al.
   "
  [term rank & {:keys [pred cnt] :or {pred identity cnt 1}}]

  (apply
   set/union
   (keep (fn[s]
           (when-let [s (filter pred (seq s))]
             (random-subset s cnt)))
         (map #(flatten
                (keep (fn[m]
                        (let [insts (taxon-direct-instances (m :taxon_id))]
                          (when insts insts)))
                      %))
              (map #(taxon-leafs (% :taxon_id))
                   (filter #(= (:node_rank %) rank)
                           (taxon-subclasses term)))))))


(defn rand-insts-by-rank-fasta-file
  "Like rand-insts-by-rank, but with the additional operation of
   taking the result set, generating an entry file from the names of
   the entries, and then generating a fasta file of the entries and
   their FULL sequences.  FILESPEC is a path to the file to be written.

   Returns full path of the fasta file generated.
   "
  [term rank filespec & {pred :pred cnt :cnt :or {pred identity cnt 1}}]

  (let [fasta-file (fs/fullpath filespec)
        entfile (fs/replace-type fasta-file ".txt")
        entries (rand-insts-by-rank term rank :pred pred :cnt cnt)]
    (edu.bc.bio.sequtils.files/gen-entry-file
     (map #(% :name) entries) entfile)
    (edu.bc.bio.sequtils.files/entry-file->fasta-file entfile)))

;;;(time (rand-insts-by-rank-fasta-file "Proteobacteria" "genus" "/data2/Bio/Training/FastaFiles/proteobacteria-genus2.fna" :pred #(= (subs (% :name) 0 2) "NC")))



;;; Main function for obtaining the bioentries requested by the user
;;; on the web site as defined by a regex they give for a taxon in the
;;; taxonomy.  Results are returned as encoded JSON, and the client
;;; issued the request as an AJAX call expecting a JSON map return
;;; value.
;;;
(defn seq-query [args]
  (let [name (args :info)
        c1 (str/take 1 name)
        name-pat (cl-format nil "[~A~A]~A"
                            c1 (str/swap-case c1) (str/drop 1 name))]
    (log> "rootLogger" :info
          "Action request: '~A', Args: '~A'"
          "Query" args)
    (let [fields "tn.name as taxname, be.*, sfcount, a.ancestors"
          tables (str "bioentry as be, taxon_name as tn, "
                      (str/replace-re #"/xxx/" name-pat
                                      taxon-nodes-table) "as tnt,"
                      seqfeature-count-table " as sfc,"
                      "ancestor as a")
          constraints ["tnt.taxon_id=tn.taxon_id"
                       "tn.name_class=\"scientific name\""
                       "tn.taxon_id=be.taxon_id"
                       "sfc.bioentry_id=be.bioentry_id"
                       "a.ncbi_taxon_id=tnt.ncbi_taxon_id"]
          constraint (str/join " and " constraints)
          stmt (str "select " fields
                    " from " tables
                    " where " constraint)]
      (sql/with-connection mysql-ds
        (sql/with-query-results qresults [stmt]
          (let [results (map xform-result
                             (doall qresults))]
            {:type :json
             :body (json/json-str
                    (conj results
                          {:me (str (str-date) ": " name-pat) :alice stmt}
                          (count results)))}
            ))))))
