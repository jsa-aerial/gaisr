;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                               A C T I O N S                              ;;
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

(ns edu.bc.bio.gaisr.actions

  "GAISR web service actions.  Compojure routes resolve to calling
   entry points here.  Generally these are 'REST'ful type actions.
   Specific result computation is performed in modules corresponding
   to the type of action requested: db querys, post processing,
   pipeline invocations, etc."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.contrib.json :as json]

            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        [edu.bc.log4clj :only [create-loggers log>]]
        [edu.bc.bio.seq-utils :only [gap-count]]
        edu.bc.bio.sequtils.files
        [edu.bc.bio.sequtils.tools :only [correct-sto-coordinates]]

        [edu.bc.bio.gaisr.db-actions
         :only [seq-query feature-query names->tax]]
        [edu.bc.bio.gaisr.pipeline
         :only [run-config-job-checked]]
        [edu.bc.bio.gaisr.post-db-csv
         :only [efile-csvhits->names-and-matches process-hit-file]]

        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]))


;;; Setup and run web based inspector.  This is pretty neat, but a bit
;;; weird as it is yet another hook into the server distinct from
;;; Emacs (or Eclipse or ...)  By default, just leave it off.
;;;
;;; NOTE: we can come in and slime-connect to any image running this
;;; and just fire it up when needed without any shutdown or restart.
;;;
;;;(use 'mycroft.main)
;;;(run 8081)


;;; Create logger(s)
;;;
;;; At the moment, only a console logger...
;;;
(create-loggers
 [[:logger "rootLogger"
   {:level :info :appenders "console"}]
  [:appender "console"
   {:type :console :pat "%-5p %d{DATE}: %m%n"}]])




;;; Support for various user actions on selected entries and items.
;;; This needs to be rethought and factored as well as it will likely
;;; become a large piece of the server functionality moving forward.
;;; Currently only supports generation of save(taxonomy, hit, ev)
;;; files, and fasta files from Scribl hit selections.
;;;
(defn save-content-xform [entries file f]
  #_(prn entries)
  (io/with-out-writer (io/file-str file)
    (doseq [e entries]
      (println e)))
  (map f entries))

(defn get-user-local-dir [user]
  (str "/home/" user "/Bio/Work/"))


(defn action-request [args]
  (let [action (args :act)
        user   (args :user)
        selections   (args :selections)]
    (log> "rootLogger" :info
          "Action Request: ~A: ~A" user action)
    (cond
     (= action "genfasta")
     (let [dir (get-user-local-dir user)]
       (cond
        (not (fs/exists? dir))
        {:type :json
         :body (json/json-str
                {:error (cl-format nil "Dir '~A' Does not exist" dir)})}

        (not (fs/writeable? dir))
        {:type :json
         :body (json/json-str
                {:error (cl-format nil "Dir '~A' Permission denied" dir)})}

        :else
        (let [fname (args :filename)
              ;; XFORM grabs the entry, but also must remove the fugly
              ;; name hack of post-db-csv/twiddle-name...
              xform #(str/replace-re #"/[0-9] " " " (subs % 0 (.indexOf % ":")))
              save-filespec (fs/fullpath (str dir fname ".txt"))
              entry-filespec (fs/fullpath (str dir fname ".ent"))
              fasta-filespec (-> (str/split #"\$\$" selections)
                                 (save-content-xform save-filespec xform)
                                 (gen-entry-file entry-filespec)
                                 (entry-file->fasta-file))]
          (log> "rootLogger" :info "~A and ~A" entry-filespec fasta-filespec)
          {:type :json
           :body (json/json-str
                  (cl-format nil "Entry File: '~A',~%Fasta File: '~A'"
                             save-filespec fasta-filespec))})))

     (= action "stocsv-match")
     (let [sto (args :sto)
           csv (args :csv)
           out (args :out)]
       (log> "rootLogger" :info " ~A<>~A --> ~A" sto csv out)
       (let [res (catch-all (efile-csvhits->names-and-matches sto csv out))]
         {:type :json
          :body (json/json-str
                 {:stat (if (list? res) "success" (str "Error: " res))})}))

     :else
     (do
       (log> "rootLogger" :info " --> ~A" selections)
       {:type :json
        :body (json/json-str
               (conj (map #(assoc {} :name %1) (str/split #" " selections))
                     {:action (str/capitalize action) :jobid 2718281828}))}
       ))))




(defn process-name2tax [filespec user-dir]
  (log> "rootLogger" :info "NAMES2TAX: ~S, ~S" filespec user-dir)
  (with-local-vars [names {}]
    (do-text-file [filespec]
      (when (and (> (count $line) 3)
                 (or (= ">gi"  (subs $line 0 3))
                     (= "#=GS" (subs $line 0 4))
                     (= "NC_"  (subs $line 0 3))))
        (var-set names (assoc (var-get names)
                         (first (re-find #"N(C|S|Z)_[0-9A-Z]+" $line)) 1))))
    (if (empty? (var-get names))
      {:error (cl-format
               nil "File '~A' Unknown format, empty, or no names"
               (fs/basename filespec))}
      (let [results (names->tax (var-get names))
            tax-filespec (fs/replace-type (fs/basename filespec) ".tax")
            tax-filespec (fs/join user-dir tax-filespec)]
        (io/with-out-writer tax-filespec
          (doseq [nv results]
            (println (str (nv :name) ": " (nv :ancestors)))))
        {:info (cl-format nil "Generated Taxonomies:~%Taxonomy File: '~A'"
                          tax-filespec)
         :tax-filespec tax-filespec}))))




(defmacro remote-try
  "Wraps a body of JSON returning remote access results in try/catches
   to ensure we don't just roll over and die with no reporting to
   client!  Accounts for the 'special' form of a contrib.condition
   Condition.
  "
  [& body]
  `(try
     (do ~@body)
     (catch clojure.contrib.condition.Condition c#
       {:body (json/json-str {:info (str (dissoc (meta c#) :stack-trace))
                              :stat "error"})})
     (catch Exception e#
       {:body (json/json-str {:info (str {:err (with-out-str (print e#))})
                              :stat "error"})})))

(defn get-remote-cemap
  "Transform the condition/exception map to a list of strings.  Remove
   stack-trace and err key/vals, place the :err value as first
   element, all others are strings of key and value as (str k \" \" v)
  "
  [cemap]
  (prn :*** cemap)
  (let [emsg (if-let [em (cemap :err)] em "Error, see below...")]
    (cons emsg
          (keep (fn[[k v]]
                  (when (not (in k [:err :stack-trace]))
                    (let [v (if (coll? v) (doall v) v)]
                      (str k " " v))))
                cemap))))


(defn remote-check-sto
  "Run check-sto on uploaded file upload-file."
  [upload-file]
  (let [result (check-sto upload-file :printem false)]
    {:body (json/json-str {:info [:NA :NA
                                  (if (= result :good)
                                   "sto is good"
                                   (cons "***BAD sto" result))]
                           :stat "success"})}))


(defn new-aln-name-seq
  "One fugly function that needs to be rather persnickety as the
   details are kind of tricky - or maybe just bean counting precise...

   Take an entry ent and its sequence sq and left/right deltas ld and
   rd and return a new [ent' sq'] pair reflecting the changes required
   by the addition/removal of characters (including gaps) from the
   front and end of sq.

   Gaps are counted as characters of sq, but not as sequence elements
   for recalculation of coordinates.
  "
  [ent sq ld rd]
  (let [sqcnt (count sq)
        [nm [os oe] sd] (entry-parts ent)
        oe (if (= oe Long/MAX_VALUE)  ; check for no coordinate info in entry
             (+ os sqcnt)
             oe)

        [ld rd sq] (if (= sd "1") ; flip for coord computation
                     [ld rd sq]
                     [rd ld (str/reverse sq)])
        s (- os (if (neg? ld) (+ ld (gap-count (str/take (- ld) sq))) ld))
        e (+ oe (if (neg? rd) (+ rd (gap-count (str/drop (+ sqcnt rd) sq))) rd))
        s (if (neg? s) 1 s) ; check for bad ld forcing falling off front
        [ld rd sq] (if (= sd "1") ; flip back for seq subs
                     [ld rd sq]
                     [rd ld (str/reverse sq)])

        zero-pos? #(>= % 0)
        prefixfn #(if (= s os)
                    ""
                    (gen-name-seq
                     (str/join "/" [nm (str s "-" (dec os)) sd])))
        suffixfn #(if (= oe e)
                    ""
                    (gen-name-seq
                     (str/join "/" [nm (str (inc oe) "-" e) sd])))

        sq (cond
            (and (zero-pos? ld) (zero-pos? rd))
            (let [[_ lssq] (prefixfn)
                  [_ rssq] (suffixfn)
                  [lssq rssq] (if (= sd "1") [lssq rssq] [rssq lssq])]
              (str lssq sq rssq))

            (and (zero-pos? ld) (neg? rd))
            (let [[_ lssq] (if (= sd "1") (prefixfn) (suffixfn))
                  rssq (str/butlast (abs rd) sq)]
              (str lssq rssq))

            (and (neg? ld) (neg? rd))
            (subs sq (abs ld) (+ sqcnt rd))

            (and (neg? ld) (zero-pos? rd))
            (let [[_ rssq] (if (= sd "1") (suffixfn) (prefixfn))
                  lssq (subs sq (abs ld))]
              (str lssq rssq)))]

    [(str/join "/" [nm (str s "-" e) sd]) sq]))


(defn remote-gen-name-seqs
  "Create a new or adjust an old set of entry/sequence pairs.  There
   are two main types of input for this:

   * Sequence files
   * Entry files

   Sequence files in turn are of two types: simple fasta format or
   alignment formats of sto or clustalw.

   Entry files are just that - simple text files of entries, one per
   line, with nc name, coordinates, and strand information. as defined
   by gen-name-seq (for which see for the details).  Entry files
   return the entries and corresponding sequences in fasta format.

   For fasta files, just treat them as an alternate form of entry file
   and return a fasta file with the adjusted sequences and
   coordinates (see delta info below).

   For alignment files, must carefully adjust the sequences and
   coordinates to account for gaps.  Deltas affect total characters
   added/removed from prefix/suffix of sequences.  However,
   coordinates for this are adjuste by only the amount of non gap
   characters in those substrings of the sequence.

   deltas is an array (possibly nil) of one or two integers (+/-) that
   indicate changes to the prefix and suffix of sequences and the
   coordinates of their corresponding entries.  Again, details are as
   defined by gen-name-seq.
  "
  [user filename upload-file deltas]
  (remote-try
   (let [alns ["sto" "aln"]
         fnas ["fna" "fa"]
         sto-n-orig-line (take 2 (io/read-lines upload-file))
         ftype (fs/ftype filename)
         [ldelta rdelta] deltas
         ldelta (if ldelta ldelta 0)
         rdelta (if rdelta rdelta 0)

         new-ents-seqs
         (if (in ftype ["ent" "txt"])
           (gen-name-seq-pairs
            (io/read-lines upload-file) :ldelta ldelta :rdelta rdelta)

           ;; Else, some sort of sequence file
           (let [entries-and-seqs (read-seqs upload-file :info :both)
                 fna (in ftype fnas)]
             (map (fn[[ent sq]]
                    (let [[nent nsq :as both]
                          (if fna ; NON alignment easy!
                            (gen-name-seq ent :ldelta ldelta :rdelta rdelta)
                            ;; Else must carefully adjust within alignment!
                            (new-aln-name-seq ent sq ldelta rdelta))]
                      both))
                  entries-and-seqs)))

         local-file (fs/tempfile "rem-gen-name-seqs-" (str "." ftype))]

     (io/with-out-writer local-file
       (when (in ftype alns)
         (let [lines sto-n-orig-line
               lines (if (= ftype "sto") (concat lines [""])  lines)]
           (doseq [l lines] (println l))))
       (doseq [[idi sq] new-ents-seqs]
         (if (in ftype alns)
           (cl-format true "~A~40T~A~%" idi sq)
           (do (println (str ">" idi)) (println sq))))
       (when (= ftype "sto")
         (doseq [gcline (filter #(.startsWith % "#")
                                (first (sto-GC-and-seq-lines upload-file)))]
           (let [[gc kind v] (str/split #" +" 3 gcline)]
             (cl-format true "~A~40T~A~%" (str gc " " kind) v)))
         (println "//")))

     {:body (json/json-str
             {:info [:NA :NA (slurp local-file)]
              :stat "success"})})))


(defn remote-correct-sto-coordinates
  "Run coordinate corrector for user on uploaded file upload-file."
  [user filename upload-file]
  (remote-try
   (let [resultfn
          (fn[[file res]]
            (if (map? res)
              (cons :error (cons file (get-remote-cemap res)))
              (cons :good (map #(if (fs/exists? %) (slurp %) []) res))))
          jobid (add *jobs* user
                     `[(print-str ~filename)
                       (correct-sto-coordinates ~upload-file)]
                     resultfn)]
      (start *jobs* jobid)
      {:body (json/json-str
              {:info [:started :started
                      (str filename " coordinate correction started, jobid ")
                      (str jobid)]
               :stat "success"})})))


(defn remote-run-config
  "Run a config job based on information in uploaded config-file"
  [user filename config-file]
  (remote-try
    (let [resultfn
          (fn[[file res]] ; input vec of task results
            (prn :### file res)
            (cond
             (map? res) (cons :error (cons file (get-remote-cemap res)))
             (not= :good res) (cons :error (cons file (seq (ensure-vec res))))
             :else (cons :good (list (str "Completed job defined by: " file)
                                     (str "Result: " res)))))
          jobid (add *jobs* user
                     `[(print-str ~filename)
                       (run-config-job-checked ~config-file :printem false)]
                     resultfn)]
      (start *jobs* jobid)
      {:body (json/json-str {:info [:started :started
                                    (str filename " job started, jobid ")
                                    (str jobid)]
                             :stat "success"})})))

(defn remote-check-job
  "Check status, and if done, return results for job JOBID of user USER"
  [user jobid]
  ;; This 'works', but it is fugly and I don't like it.  Needs to be
  ;; refactored and cleaned up.
  (remote-try
   (let [stat (catch-all (check-job jobid))
         res (if (= stat :done)
               (seq (ensure-vec (result *jobs* jobid)))
               (str "currently "
                    (if (= :in-process (result *jobs* jobid :t2))
                      "running"
                      (str (result *jobs* jobid :t2)))))
         [jstat res] (if (= stat :done)
                       [(first res) (rest res)]
                       [:running res])]
     {:body (json/json-str
             {:info (concat
                     [stat jstat (str "Job " jobid " for user '" user "'")]
                     (ensure-vec res))
              :stat "success"})})))


(defn upload-file
  "Dispatch function for upload posts.  ARGS contains a map of
   subcommand information, including the primary function to perform
   encoded as value of key upload-type.
  "
  [args reqmap]
  (log> "rootLogger" :info "UPLOAD ~S, Cookies: ~S" args (reqmap :cookies))
  (let [user (get-in reqmap [:cookies "user" :value])
        upload-type (args :upload-type)
        fileinfo (args :file)
        upload-file (fileinfo :tempfile)
        filename (fileinfo :filename)
        local-out (fs/fullpath (str "/data2/Bio/Work/" filename))]
    (cond
     (= filename "")
     {:body (json/json-str {:error "No file specified"})}

     (= (fileinfo :size) 0)
     {:body (json/json-str
             {:error (cl-format nil "File '~A' has no content" filename)})}

     :else
     (cond
      (= upload-type "hitfile")
      (do (io/copy upload-file (io/file-str local-out))
          {:body (json/json-str (process-hit-file local-out))})

      (= upload-type "check-sto")
      (remote-check-sto upload-file)

      (= upload-type "correct-sto-coordinates")
      (let [tmp-sto (fs/replace-type (.getPath upload-file) ".sto")]
        (io/with-out-writer tmp-sto (print (slurp upload-file)))
        (remote-correct-sto-coordinates user filename tmp-sto))

      (= upload-type "get-seqs")
      (let [deltas (map #(if (string? %) (Integer. %) %) (args :misc))
            ftype (fs/ftype filename)
            tmp-file (fs/replace-type (.getPath upload-file) (str "." ftype))]
        (fs/rename (.getPath upload-file) tmp-file)
        (remote-gen-name-seqs user filename tmp-file deltas))

      (= upload-type "run-config")
      (remote-run-config user filename (.getPath upload-file))

      (= upload-type "check-job")
      (remote-check-job user (read-string (slurp upload-file)))

      (= upload-type "name2tax")
      {:body (json/json-str
              (cond
               (nil? user)
               {:error "File to Taxonomy file requires USER"}

               (= (args :subtype) "remote")
               (do (io/copy upload-file (io/file-str local-out))
                   (let [udir (get-user-local-dir user)
                         txf (:tax-filespec (process-name2tax local-out udir))]
                     {:info (slurp txf) :stat "success"}))
               :else
                (do (io/copy upload-file (io/file-str local-out))
                    (process-name2tax local-out (get-user-local-dir user)))))}

      (= upload-type "stocsv-match")
      (let [csvinfo (args :csv)
            csv-upload (csvinfo :tempfile)
            csv-filename (csvinfo :filename)
            sto local-out
            csv (fs/fullpath (str "/data2/Bio/Work/" csv-filename))
            tmp (fs/tempfile "stocsvout" ".tmp" "/data2/Bio/Work/")]
        (io/copy upload-file (io/file-str sto))
        (io/copy csv-upload (io/file-str csv))
        (let [res (catch-all (efile-csvhits->names-and-matches sto csv tmp))]
          {:type :json
           :body (json/json-str
                  (if (list? res)
                    {:info (slurp tmp) :stat "success"}
                    {:stat (str "Error: " res)}))}))

      :else
      {:body (json/json-str
              {:error "'Dunno yet' means don't know yet! :-)"})}))))




;;; Wrapper for module/name-space separation and visibility structure
;;;
(defn get-features [args]
  (feature-query args))


;;; Wrapper for module/name-space separation and visibility structure
;;; Obsolete?  Check and remove!
(defn map-names-to-ancestors [args]
  (names->tax args))


;;; Sequence query by context.  Currently just by taxon - does not yet
;;; account for the new expansion in querying capabilities
;;; incorporating genomic context.  This will support queries taking
;;; into account constraints of seqfeature types and values (eg.,
;;; gene=rplK), locus (uploc=5000 and downloc=3000), et.al.
;;;
(defn ctx-query [args reqmap]
  (let [response (seq-query (merge {:info (args :query)} args))]
    (if (get-in reqmap [:cookies "user" :value])
      response
      (assoc response
        :cookies {:user {:value (:user args)
                         :max-age (str (* 3600 24 30 12))}}))))




;;; Default request handler for non existent routes
;;;
(def req-body (atom ""))

(defn fall-through-catcher [uri params form-params request]
  (log> "rootLogger" :warn
        "URI '~A'~%Method: ~A~%Params: ~A~%Form-parms: ~A~%~%~A"
        uri (request :request-method) params form-params request)
  (when (request :body)
    (swap! req-body (fn[_ rb] rb) (request :body)))
  (str "<h1>Unknown Request: '" uri
       "', with query info '" (request :query-string)
       "'</h1>"))
