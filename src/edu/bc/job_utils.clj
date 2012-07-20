(in-ns 'edu.bc.utils)


;;; ----------------------------------------------------------------
;;; Long term job/thread/task run management with futures.  This may
;;; be subsumable by a "reasonable" batch system (beanstalked wrapped
;;; via Sisyphus or something).  However, most of such stuff looks
;;; like overkill, but maybe that will change.  If so, perhaps revisit
;;; and use the "reasonable" thing...


(defrecord job
  [id tasks status results])


(defn make-job [jid task-set]
  (job. jid (atom task-set) (atom {jid :new}) (atom {jid :none})))


;;; Too cute, but worse - no cycle check.
;;; (defn job-top-sortx [x]
;;;   (mapcat #(for [[k v] %1 :when (empty? v)] k)
;;;       (take-while seq (iterate
;;;                        #(into {} (for [[k v] % :when (seq v)]
;;;                                    [k (mapcat % v)]))
;;;                        x))))

(defprotocol job-stuff
  "Protocols defining the the character of a jobs maintenance and processing
   capability.  The main functions are adding new jobs, starting them, and
   interogating them for results (partial as well as full) and on going status
   of tasks as they run and complete, and the final job results.

   Job addition builds an implicit constraint graph where each node is a
   task of the job (a separate function call where some of the arguments
   may be other task results).  Each such specification is transformed into
   a call against the arguments, where any task argument is transformed into
   a dereference of the task's result.  This call C is then used as the main
   processing of a computed function for the task which creates a promise
   in the job's results map indexed by the task's id and a future which will
   set start status, fire a future wrapping the body which delivers the
   result of the body to the promise.  When called this function fires the
   future and returns it.  Since the arguments that are task ids in call
   C turn to references of the tasks result promises, a dependency data
   flow graph is created.  Tasks run in parallel unless they require the
   results of other task, in which case they wait until those results are
   available.  The graph is a DAG and enforced as such."

  (all
   [jobs]
   [jobs user]
   [jobs user stat]
   "No parameters: Return a seqable of all the current job ids.
    With STAT, return a seqable of all job ids with status STAT.
    User parameter: Return a seqable of all the known jobs for USER.
    If STAT return all known jobs for USER with status STAT.
    In all cases, this may include those jobs saved in the database.")

  (all-tasks [jobs jid] "Return all task ids for job with id jid")

  (add
   [jobs user tasks jobfn]
   "Add new job for user USER with task set tasks and final job 'reduction'
    jobfn.  Add new job's id to those for USER.  Return new job's ID

    TASKS is a seq (typically a vector) of forms (typically function calls)
    which perform the tasks of the job.  Could be only one, but there must be
    at least one.  NOTE: these are forms, i.e. they must be unevaluated and
    thus must be quoted.

    Tasks can depend on the results of other tasks, but there cannot be a
    cycle (the dependency graph is a DAG).  A topological sort is performed
    at task set creation time and any cycle will raise and exception.

    A dependency is specified by using a keyword :tn, n an integer in
    1..(count tasks), as an argument in the form of a task.  For
    example, if the task form is a function call, one or more of its
    arguments can refer to some other task result as :tn.  At task
    construction time, these argumets are transformed into promises for
    the specified task values with the result being a dataflow constraint
    network.  All tasks whose arguments are available will run in parallel.

    JOBFN is the 'reduction' function for the job.  It must be a single
    parameter function and it is passed a seq of the results of the tasks
    (in order).  It's result is the result value of the job.
    ")

  (fetch [jobs jid] "Fetch the job record with id jid")
  (start [jobs jid] "Start a job.  Launches all task futures")

  (chk
   [jobs jid]
   [jobs jid tid]
   "check on status of job JID or task TID of job with JID")

  (result
   [jobs jid]
   [jobs jid tid]
   "Return result of job JID or task TID of job JID")

  (wresult [jobs jid] "Wait for, and then return, the result of job JID")
  (results [jobs jid] "Return the results of all tasks and the job")
  (cancel [jobs jid] "cancel job with id JID.")

  (del [jobs user jid & {:keys [force] :or {force false}}]
       "remove job with id JID from job set"))




;;; Job task graph topological sort of dependencies.
;;;
(defn next-task [deps done]
  (some (fn [[k [_ ts] :as task]]
          (when (empty? (remove done ts)) task))
        deps))

(defn job-top-sort [tasks]
  (loop [deps tasks
         done #{}
         tlist []]
    (if (empty? deps)
      tlist
      (if-let [[k [fn _]] (next-task deps done)]
        (recur (dissoc deps k)
               (conj done k)
               (conj tlist [k fn]))
        (raise :type :circular-task-dependency
               :tasks tasks
               :tlist tlist
               :state [deps done])))))


;;; Used as runtime dereference of tasks result promises of dependent
;;; tasks
;;;
(defn task-result [j tid]
  ;;(println "TaskRes" tid (deref (tid @(:results j))))
  (deref (tid @(:results j))))


(declare *jobs*)

;;; Main task set constructor.  Builds and compiles a function that is
;;; the body of each task taking into account status and result
;;; updates, and proper waiting on required upstream results.  Also
;;; builds and compiles the body of the job reduction function where
;;; it waits on all task results.
;;;
(defn- make-task-functions [tasks jobfn]
  (let [jid (gen-kwuid)
        tcnt (count tasks)
        taskids (map #(keyword (str %1 %2)) (repeat tcnt "t") (iterate inc 1))

        ;; Generate the dependency map for the job topological sort.
        ;; The entries are: tid [task-func dependencies].
        dep-map
        (apply merge
               (map
                (fn [task tid]
                  (let [task-args (keep #(when (in %1 taskids) %1)
                                        (rest task))

                        ;; Here we filter the arguments so that any
                        ;; which refer to a task (i.e, a keyword of
                        ;; the form :ti) are changed to reference the
                        ;; value of the promise of the corresponding
                        ;; task
                        deref-args (map #(if (in %1 taskids)
                                           `(task-result ~'j ~%1)
                                           %1)
                                        (rest task))

                        ;; Now build the new call which will be the
                        ;; task function body.
                        body `(~(first task) ~@deref-args)

                        ;; Make and compile a function fn for this
                        ;; task. fn creates a promise and associates
                        ;; it with the task id in the results of our
                        ;; job.  fn then fires a future which delivers
                        ;; the result of the form passed as the task
                        ;; body to this promise.  fn returns the
                        ;; _future_.
                        task (eval
                              `(fn[]
                                 (let [~'p (promise)
                                       ~'j (fetch *jobs* ~jid)]
                                   ;;(println "Starting:" ~tid)
                                   (swap! (:results ~'j) assoc ~tid ~'p)
                                   (future
                                    (swap! (:status ~'j) assoc ~tid :waiting)
                                    (let [~'args [~@deref-args]]
                                      (swap! (:status ~'j) assoc ~tid :running)
                                      (let [~'v (deliver
                                                 ~'p (catch-all ~body))]
                                        (swap! (:status ~'j)
                                               assoc ~tid :done)
                                        ~'v))))))]
                    {tid [task task-args]}))
                tasks taskids))

        ;; Special job "task" is "reduction/combination" of task
        ;; results by means of user passed in function applied to all
        ;; the task results
        jfn (eval `(fn[]
                     (let [~'p (promise)
                           ~'j (fetch *jobs* ~jid)]
                       (swap! (:results ~'j) assoc ~jid ~'p)
                       (future
                        (swap! (:status ~'j) assoc ~jid :waiting)
                        (let [~'args [~@(map #(do `(task-result ~'j ~%1))
                                             taskids)]]
                          (swap! (:status ~'j) assoc ~jid :running)
                          (let [~'v (deliver
                                     ~'p (~jobfn ~'args))]
                            (swap! (:status ~'j) assoc ~jid :done)
                            ~'v))))))

        ;; Include jfn in dependencies map
        dep-map (assoc dep-map jid [jfn taskids])

        ;; Linearize the constraint network.  The reasons we need this
        ;; are 1. to raise an exception if there is a cycle and
        ;; 2. ensure the _startup_ order works wrt promises refed in
        ;; the body function call of the futures.  The promise
        ;; constraint graph will take care of the _processing_
        ;; ordering automagically.
        task-seq (job-top-sort dep-map)]

    [jid task-seq]))




(defn make-task-runner
  "Takes TASK-SEQ which is a set of [tid tfn] pairs and builds a function
   which maps the pairs to [tid (tfn)], i.e., the ids paired with the results
   of calling the tfns, and inserts the result into a map and returns it.
   Each tfn launches and returns a future.  The resulting map is used as
   the corresponding job's tasks."
  [task-seq]
  (fn [_]
    (into {} (map (fn[[tid tfn]] [tid (tfn)]) task-seq))))


(defmacro with-base-chk [jid & body]
  `(let [jid# ~jid
         j# (@~'jmap jid#)]
     (if (nil? j#)
       (raise :type :no-such-job :id jid#
              :msg (cl-format nil "No active job with id ~A" jid#))
       (do ~@body))))


(defrecord jobs [jmap umap]
  job-stuff
  (all [_] (map (fn[[_ job]] (:id job)) @jmap))
  (all [_ user]
       (@umap user []))
  (all [this user stat]
       (let [ujobs (@umap user [])]
         (keep #(when (= stat (chk this %1)) %1) ujobs)))

  (all-tasks [jobs jid] (map (fn[[tid _]] tid) @(:tasks (fetch jobs jid))))

  (add [_ user tasks jobfn]
       (let [[jid task-seq] (make-task-functions tasks jobfn)
             trun (make-task-runner task-seq)
             j (make-job jid {:trun trun})]
         (swap! jmap assoc jid j)
         (swap! umap (fn[m u jid]
                       (assoc m u (conj (get m u []) jid)))
                user jid)
         jid))

  (fetch [_ jid]
         (@jmap jid))

  (start [_ jid]
         (with-base-chk jid
           (let [j (@jmap jid)
                 trun (@(:tasks j) :trun)]
             (assert (fn? trun))
             (swap! (:tasks j) trun))))

  (chk [_ jid]
       (with-base-chk jid
         (@(:status (@jmap jid)) jid)))
  (chk [_ jid tid]
       (with-base-chk jid
         (@(:status (@jmap jid)) tid :not-available)))

  (result [_ jid]
          (with-base-chk jid
            (case (chk _ jid)
                  :waiting :waiting
                  :running :in-process
                  :new :not-started
                  :done @(@(:results (@jmap jid)) jid))))
  (result [_ jid tid]
          (if (= jid tid)
            (result _ jid)
            (with-base-chk jid
              (case (chk _ jid tid)
                    :waiting :waiting
                    :running :in-process
                    :not-available :not-started
                    :done @(@(:results (@jmap jid)) tid)))))

  (wresult [this jid]
           (with-base-chk jid
             (let [j (fetch this jid)]
               ;; First, make sure the job is started!
               (while (= :new (chk this jid))
                 (Thread/sleep 200))
               (task-result j jid))))

  (results [this jid]
           (with-base-chk jid
             (into {} (map (fn[[k _]]
                             [k (if (= k jid)
                                  (result this k)
                                  (result this jid k))])
                           @(:results (@jmap jid))))))

  (cancel [_ jid]
          (with-base-chk jid
            (raise :type :nyi :msg "job cancel not yet implemented!")))

  (del [_ user jid & {:keys [force] :or {force false}}]
       (with-base-chk jid
         (when (not force)
           (assert (in (chk _ jid) [:new :waiting :done :cancelled :abort])))
         (dosync
          (swap! jmap dissoc jid)
          (swap! umap (fn[m u jid]
                        (assoc m u (remove  #(= %1 jid) (get m u))))
                 user jid)))))


;;; Map of all jobs this server run (i.e., non persistent, but
;;; complete).  It isn't clear that keeping information about old done
;;; jobs is useful and thus makes any sense.  If it does, then really
;;; this should likely be moved to the database (typical transactional
;;; schema).  But until it looks like that makes sense, just use a map
;;; of job ids to run set, promises for tracking results, and futures
;;; moving along the ongoing processes of the job.
;;;
;;; A "job" here is just any request for "batch" type processing from
;;; someone.  And "batch" means something loosely like "takes a while
;;; and/or has several separate parts (typically with dependencies
;;; between them) to the overall processing".
;;;

(def *jobs* (jobs. (atom {}) (atom {})))



(defn all-jobs
  ([stat]
     (keep #(when (= stat (chk *jobs* %1)) %1)
           (all *jobs*)))
  ([user stat]
     (if (= stat :all)
       (all *jobs* user)
       (all *jobs* user stat))))


(defn- start-set [jids]
  (map #(start *jobs* %1) jids))

(defn start-all
  ([] (start-set (all *jobs* :new)))
  ([user]
     (start-set (all *jobs* user :new))))


(defn- delete-set [user jids force]
  (map #(del *jobs* user %1 :force force) jids))

(defn delete-done [user]
  (delete-set user (all *jobs* user :done) false))

(defn delete-all [user & {:keys [force] :or {force false}}]
  (delete-set user (all *jobs* user) force))


(defn- check-set [jids]
  (map #(chk *jobs* %1) jids))

(defn check-jobs
  ([] (check-set (all *jobs*)))
  ([user]
     (check-set (all *jobs* user))))

(defn check-job
  [jobid]
  (chk *jobs* jobid))

(defn check-tasks [jid]
  (let [jids (if (coll? jid) jid [jid])]
    (map (fn[jid]
           (map #(chk *jobs* jid %1)
                (all-tasks *jobs* jid)))
         jids)))

(defn check-all
  ([]
     (map #(check-tasks %1)
          (sort (all *jobs*))))
  ([user]
     (map #(check-tasks %1)
          (sort (all *jobs* user)))))


(defn- results-for-set [s kind]
  (if (= kind :jobs)
    (map #(do [%1 (result *jobs* %1)]) s)
    (map #(results *jobs* %1) s)))

(defn results-jobs
  ([]
     (results-for-set (sort (all *jobs*)) :jobs))
  ([user]
     (results-for-set (sort (all *jobs* user)) :jobs)))

(defn results-all
  ([]
     (results-for-set (sort (all *jobs*)) :all))
  ([user & stat]
     (cond
      (nil? stat)
      (results-for-set (sort (all *jobs* user)) :all)

      (not (coll? stat))
      (results-for-set (sort (all *jobs* user stat)) :all)

      :else
      (let [jids (sort (flatten
                        (concat (map #(all *jobs* user %1)
                                     stat))))]
        (results-for-set jids :all)))))

