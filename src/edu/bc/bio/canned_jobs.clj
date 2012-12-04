;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                           C A N N E D - J O B S                          ;;
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

(ns edu.bc.bio.canned-jobs

  "General utility functions and macros.  Basically these resources are
   fairly general and intended to be usable on most any project"

  (:require [clojure.contrib.math :as math]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.tools

        [edu.bc.log4clj :only [create-loggers log>]]
        [clojure.pprint
         :only [cl-format]]))



;;; ------------------------------------------------------------------------;;;
;;; DB loading jobs and job support

(defn gunzip-file [dir file]
  (let [dir (if (= (last dir) \/) dir (str dir "/"))]
    (catch-all (runx "gunzip" (str dir file)))))

(defn load-gbff-file [filespec]
  (let [filein (fs/fullpath filespec)
        bp-path (get-tool-path :bioperl)
        bp-load-seqdb (str bp-path "bp_load_seqdatabase.pl")
        dbconnect-args ["--driver" "mysql" "--host" "127.0.0.1" "--port" "3306"]
        dbuser-args ["--dbuser" "root" "--dbpass" "rna314rulz"]
        dbname ["--dbname" "biosql"]
        ctrl-args ["--format" "genbank" "--safe" "--noupdate" "--lookup"]
        switches (concat dbconnect-args dbuser-args dbname ctrl-args)
        cmd-args (conj switches filein)]
    (assert-tools-exist [bp-load-seqdb])
    (catch-all (runx bp-load-seqdb cmd-args))))


;;; Tasks 1 and 2 for genbank file format files
;;;
(defn t1-gbff-load [dir file]
  (let [result (gunzip-file dir file)]
    (if (string? result)
      (str dir (str/replace-re #".gz" "" file))
      [:failed result])))

(defn t2-gbff-load [filespec & _]
  (if (not (string? filespec))
    [:failed filespec]
    (let [result (catch-all (load-gbff-file filespec))]
      (if (string? result)
        (str filespec " loaded")
        [:failed result]))))


;;; For now we are still using bioperl db loader.  This is not
;;; thread/process safe and so we must force serialization of loader
;;; jobs as well as tasks in them.  SERIALIZE-LOAD will be a task in
;;; each job (except the first!) in a set of jobs that load gbff files
;;; to the db.  This will force each job to wait on it's predecssor -
;;; the scheduler sets up the job ids in order and hands them out so
;;; each successor knows its predecessor.  When we move to a Clojure
;;; loader, we will be able to dispense with this and run any number
;;; of loaders in parallel!
;;;
(defn serialize-load [jid]
  (wresult *jobs* jid))


;;; The job result function for a load.  ARGS is the set of results of
;;; all the tasks for the job in topological order.  Just report
;;; finished if they all succeeded, or failed otherwise for the job's
;;; result.
;;;
(defn jfn-load [args]
  (let [res (reduce #(if (and (string? %1) (string? %2) (not= %1 "FAILED"))
                       "FINISHED"
                       "FAILED")
                    "" args)]
    (println (str-date) (if (> (count args) 2) (rest args) args))
    res))


;;; Genbank file format file loader scheduler.  Take CNT gbff files
;;; from directory, that have yet to be loaded (defined as still
;;; gzipped) and create a task schedule for each which will load the
;;; file and for all but the first, also add a serializer task so that
;;; the scheduled jobs will run serially.  Needed due to current use
;;; of non thread/process safe bioperl loader...
;;;
(defn schedule-gb-files [directory cnt]
  (let [directory (fs/fullpath directory)
        directory (if (= \/ (last directory)) directory (str directory "/"))
        dircontent (catch-all (runx "ls" directory))]
    (when (map? dircontent)
      (raise dircontent))
    (when (not (string? dircontent))
      (raise :type :unknown-result
             :cmd (str "ls " directory)
             :result dircontent))
    (let [file-names (take cnt (sort (keep #(when (re-find #"\.gz$" %1) %1)
                                           (str/split #"\n" dircontent))))]
      (loop [fnames file-names
             id nil]
        (if (empty? fnames)
          (conj file-names :scheduled)
          (let [fname (first fnames)
                jid (if id
                      (add *jobs*
                           "mlab"
                           `[(serialize-load ~id)
                             (t1-gbff-load ~directory ~fname)
                             ; wait for file & prev load.  The opacity
                             ; of this is probably a good indication
                             ; to not use job dependencies for current
                             ; non thread safe bioperl loader!
                             (t2-gbff-load :t2 :t1)]
                           'jfn-load)
                      (add *jobs*
                           "mlab"
                           `[(t1-gbff-load ~directory ~fname)
                             (t2-gbff-load :t1)]
                           'jfn-load))]
            (recur (rest fnames) jid)))))))




;;; ------------------------------------------------------------------------;;;
;;; Pipeline jobs.

(def default-work-dir
     (or (getenv "BIO_DEFAULT_WORK_DIR")
         "~Bio/Work/"))

(defn get-work-dir [& [dir]]
  (let [dir (fs/fullpath (if dir dir default-work-dir))
        dir (if (= \/ (last dir)) dir (str dir "/"))]
    dir))


(defn get-range-stg [range]
  (->> range
       (#(if (string? %) % ""))
       (str/replace-re #"/-1" "-minus")
       (str/replace-re #"/1" "-plus")))

(defn get-entry-file-name [entries]
  (let [[entry range] (str/split #" " (first entries))
        ns (str entry (get-range-stg range))
        [entry range] (str/split #" " (last entries))
        ne (str entry (get-range-stg range))
        fname (str ns "-" ne ".txt")]
    fname))



;;; Blasting: From entry file generation to parsing out for "relevant"
;;; matches.  For mlab, typically all pieces of this will be used for
;;; a given run.
;;;
(defn t-gen-entry-file [entries outfile]
  (catch-all
   (gen-entry-file entries (fs/fullpath outfile))))

(defn t-entry-file->fasta-file [efile]
  (catch-all
   (entry-file->fasta-file efile)))


(defn t-blastn [in args]
  (catch-all (apply blastn in args)))

(defn t-tblastn [in args]
  (catch-all (apply tblastn in args)))


(defn t-blast-parse [in & [out]]
  (catch-all
   (let [in (fs/fullpath in)
         out (if out (fs/fullpath out) (str in ".matches"))]
     (parse-blast in out))))


(defn jfn-blast [args]
  (map #(catch-all (runx "ls" "-l" %1)) args)
  (last args))


(defn schedule-blast-seq [user pgm entries dir args]
  (let [dir (get-work-dir dir)
        out-filespec (str dir (get-entry-file-name entries))]
    (add *jobs*
         user
         `[(t-gen-entry-file '~entries ~out-filespec)
           (t-entry-file->fasta-file :t1)
           (~pgm :t2 ~@args)
           (t-blast-parse :t3)]
         'jfn-blast)))

(defn schedule-blast-ent [user pgm efile dir args]
  (let [efile (fs/fullpath efile)]
    (add *jobs*
         user
         `[(t-entry-file->fasta-file ~efile)
           (~pgm :t1 ~@args)
           (t-blast-parse :t2)]
         'jfn-blast)))

(defn schedule-blast-fna [user pgm fnafile dir args]
  (add *jobs*
       user
       `[(~pgm ~fnafile ~@args)
         (t-blast-parse :t1)]
       'jfn-blast))

(defn schedule-blast
  [user entries &
   {pgm :pgm etype :etype dir :dir args :args
    :or {pgm (blastpgm (get-selection-fna entries))
         etype :fna dir (get-work-dir)
         args [:word-size (if (= pgm tblastn) 4 8)]}}]
  (let [bfn (case etype
                  :seq schedule-blast-seq
                  :ent schedule-blast-ent
                  :fna schedule-blast-fna)]
    (bfn user pgm entries dir args)))


(defn start-scheduled-blast [jid]
  (start *jobs* jid))

(defn run-blast [entries & [dir]]
  (let [jid (schedule-blast entries dir)]
    (start *jobs* jid)
    jid))




;;; Generating a binary blastdb from a set of sequence entries.
;;;
(defn gen-blastdb [entries dir]
  (let [dir (get-work-dir dir)
        efile (gen-entry-file entries (str dir (get-entry-file-name entries)))]
    (entry-file->blastdb efile)))

(defnk run-gen-blastdb
  "Run a 'generate binary blastdb' job.  USER is passwd file username,
   ENTRIES is a coll of strings of the form 'GNAME [s-e[/(-1|1)]]'
   where GNAME is the name of a genome, s-e is an _optional_ range,
   where S = start index and E = end index for the genome sequence,
   and /(-1|1) is an _optional_ strand designator, where -1 is for
   minus and 1 is for plus.  Default range is whole sequence.  Default
   strand is 1 (plus).

   DIR is an optional working director for the output, defaults to ~Bio/Work.
   BLASTDB is an optional NCBI binary entry DB, defaults to other_genome.

   Example entries: NC_000964, NC_000964 119-777, NC_000964 119-777/-1"
  [user entries :dir (get-work-dir) :blastdb nil]
  (let [loc (some #(re-find #"[0-9]+\-[0-9]+" %) entries)
        jid (add *jobs*
                 user
                 `[(gen-blastdb ~entries ~dir blastdb loc)]
                 #(first %1))]
    (start *jobs* jid)))


;;; ---------------------------------------------------
