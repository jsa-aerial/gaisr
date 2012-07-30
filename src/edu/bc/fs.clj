(ns edu.bc.fs
  ^{:doc "File system utilities in Clojure"
    :author "Miki Tebeka <miki.tebeka@gmail.com>"}
  (:refer-clojure :exclude [empty?])
  (:require [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [clojure.zip :as zip])
  (:import java.io.File
           java.io.FileInputStream
           java.io.FileOutputStream
           java.io.FilenameFilter))

;;; (compile 'edu.bc.fs)


(def separator File/separator)


(defn fullpath
  "Canonicalize FILESPEC, a string, to a fully qualified file path
   for the native system.  ~ in position one is translated to the
   users home directory path, / and \\ are translated to the file
   separator for the native system."
  [filespec-str]
  (let [^String s (str filespec-str)
        s (.replace s  \\ File/separatorChar)
        s (.replace s \/ File/separatorChar)]
    (if (.startsWith s "~")
      (str (System/getProperty "user.home")
           separator
           (subs s 1))
      s)))


(declare exists? size)

(defn empty? [path]
  "Returns false if either (io/file path) does not exist OR if the
   denoted file is empty (has size 0)"
  (or (not (exists? path))
      (= (size path) 0)))


(defn pwd []
  (io/pwd))


(defn homedir []
  (System/getProperty "user.home"))




(defn join
  "Join parts of path.\n\t(join [\"a\" \"b\"]) -> \"a/b\""
  [& parts]
  (apply str (interpose separator parts)))

(defn split
  "Split path to componenets.\n\t(split \"a/b/c\") -> (\"a\" \"b\" \"c\")"
  [path]
  (into [] (.split path separator)))

(defn rename
  "Rename old-path to new-path."
  [old-path new-path]
  (.renameTo (io/file old-path) (io/file new-path)))




(defn exists?
  "Return true if path exists."
  [path]
  (.exists (io/file path)))

(defn directory?
  "Return true if path is a directory."
  [path]
  (.isDirectory (io/file path)))

(defn file?
  "Return true if path is a file."
  [path]
  (.isFile (io/file path)))

(defn executable?
  "Return true if path is executable."
  [path]
  (.canExecute (io/file path)))

(defn readable?
  "Return true if path is readable."
  [path]
  (.canRead (io/file path)))

(defn writeable?
  "Return true if path is writeable."
  [path]
  (.canWrite (io/file path)))




;;; From clj-file-utils
(defn delete
  "Delete path."
  [path]
  (.delete (io/file path)))

(defn rm
  "Remove a file. Will throw an exception if the file cannot be deleted."
  [file]
  (io/delete-file file))

(defn rm-f
  "Remove a file, ignoring any errors."
  [file]
  (io/delete-file file true))

(defn rm-r
  "Remove a directory. The directory must be empty; will throw an exception
if it is not or if the file cannot be deleted."
  [path]
  (io/delete-file-recursively path))

(defn rm-rf
  "Remove a directory, ignoring any errors."
  [path]
  (io/delete-file-recursively path true))




(defn abspath
  "Return absolute path."
  [path]
  (.getAbsolutePath (io/file path)))

(defn- strinfify [file]
  (.getCanonicalPath file))

(defn normpath
  "Return nomralized (canonical) path."
  [path]
  (strinfify (io/file path)))

(defn basename
  "Return basename (file part) of path.\n\t(basename \"/a/b/c\") -> \"c\""
  [path]
  (.getName (io/file path)))

(defn dirname
  "Return directory name of path.\n\t(dirname \"a/b/c\") -> \"/a/b\""
  [path]
  (.getParent (io/file path)))

(defn ftype
  "Return the file type suffix of PATH.  File type suffix is the last
   substring element following a '.' in (basename PATH).  If no such
   element is in the basename, returns empty string\"\"

   Ex: (ftype \"/home/fred/Bio/test.sto\") ==> \"sto\"
  "
  [path]
  (let [bn (basename path)
        ft (last (str/split #"\." bn))]
    (if (= ft bn) "" ft)))

(defn replace-type
  "Replace the file extension type of FILESPEC to be EXT.  The type
   for FILESPEC is the last part dotted extension.  Formally, matches
   regexp '\\.[^.]*$'.  If EXT is a seq/vec replace extensions in
   last/first pairings.  Last extension replace by (first EXT), then
   last of that result is replaced by (second EXT), etc."
  [filespec ext]
  (let [rep-type (fn [filespec ext]
                   (let [dir (dirname filespec)
                         fname (str/replace-re #"\.[^.]*$" ext
                                               (basename filespec))]
                     (if dir (str dir separator fname) fname)))]
    (reduce #(rep-type %1 %2) filespec (if (coll? ext) ext [ext]))))



(defn mtime
  "Return file modification time."
  [path]
  (.lastModified (io/file path)))

(defn size
  "Return size (in bytes) if file."
  [path]
  (.length (io/file path)))




(defn listdir
  "List files under path."
  [path]
  (seq (.list (io/file path))))

(defn- _re-dir-files [directory pat]
  (map #(join directory %)
       (filter #(re-find pat %) (listdir directory))))

(defn- fix-file-regex [l]
  (-> l str/trim (str "$") ((partial str/replace-re #"\*" ".*"))))

(defn directory-files
  "Return full path qualified file specifications for all files in
   directory whose suffix matches file-type.  Typically file-type
   would be the type suffix (e.g., .txt or .sto or .fna or .clj or
   whatever), but could include any part of the file name's suffix.
   So, giving -new.sto for example, works as well.
  "
  [directory file-type]
  (_re-dir-files directory (re-pattern (str file-type "$"))))

(defn re-directory-files
  "Return full path qualified file specifications for all files in
   directo ry whose name is matched by re.  RE is a regexp (#\"regex
   def\" literal or a string that defines a regexp (which will be
   turned into a pattern).
  "
  [directory re]
  (_re-dir-files directory (re-pattern (fix-file-regex re))))




(defn mkdir
  "Create a directory."
  [path]
  (.mkdir (io/file path)))

(defn mkdirs
  "Make directory tree."
  [path]
  (.mkdirs (io/file path)))

(defn copy [from to]
  (let [from (io/file from)
        to (io/file to)]
    (when (not (.exists to)) (.createNewFile to))
    (with-open [to-channel (.getChannel (FileOutputStream. to))
                from-channel (.getChannel (FileInputStream. from))]
      (.transferFrom to-channel from-channel 0 (.size from-channel)))))


(defn dodir
  "Map ACTIONF over application of FILTERF to DIR.

   FILTERF is a function that operates on a directory returning a
   result set (typically a set of files).  ACTIONF is a function that
   takes an element of the set and any ancillary supplied arguments
   ARGS and returns an appropriate result.

   Returns the seq of results of ACTIONF minus any nils"
  [dir filterf actionf & args]
  (let [files (filterf dir)]
    (keep #(when-let [r (apply actionf % args)] r) files)))


(defn cpfiles
  "Copy a set of files from a set of dirs DIRDIR to an output dir OUTDIR.

   FILE-TYPE is the type of files in each dir in DIRDIR that are
   candidates to be copied.  REGEX is a regular expression to filter
   the set of candidates

   EX: (cpfiles \"/data2/Bio/Training/MoStos2\" \"/home/jsa/TMP/X\"
                 #\"firm\" \".sto\")
  "
  [dirdir outdir regex file-type]
  (flatten
   (dodir dirdir
          #(directory-files % "")
          (fn[d]
            (dodir d #(directory-files % file-type)
                   #(when (re-find regex %)
                      (let [to (join outdir (basename %))]
                        (copy % to) to)))))))



; FIXME: Write this
; (defn copytree [from to] ...

(defn tempfile
  "Create a temporary file."
  ([] (tempfile "-fs-" ""))
  ([prefix] (tempfile prefix ""))
  ([prefix suffix] (.getAbsolutePath (File/createTempFile prefix suffix)))
  ([prefix suffix directory]
   (.getAbsolutePath (File/createTempFile prefix suffix (File. directory)))))

(defn tempdir
  "Create a temporary directory"
  ([] (let [dir (File/createTempFile "-fs-" "")
            path (.getAbsolutePath dir)]
        (.delete dir)
        (.mkdir dir)
        path))
  ([root]
   (let [dir (File/createTempFile "-fs-" "" (File. root))
         path (.getAbsolutePath dir)]
     (.delete dir)
     (.mkdir dir)
     path)))

(defn cwd
  "Return the current working directory."
  []
  (abspath "."))




; Taken from https://github.com/jkk/clj-glob. (thanks Justin!)
(defn- glob->regex
  "Takes a glob-format string and returns a regex."
  [s]
  (loop [stream s
         re ""
         curly-depth 0]
    (let [[c j] stream]
        (cond
         (nil? c) (re-pattern (str (if (= \. (first s)) "" "(?=[^\\.])") re))
         (= c \\) (recur (nnext stream) (str re c c) curly-depth)
         (= c \/) (recur (next stream) (str re (if (= \. j) c "/(?=[^\\.])"))
                         curly-depth)
         (= c \*) (recur (next stream) (str re "[^/]*") curly-depth)
         (= c \?) (recur (next stream) (str re "[^/]") curly-depth)
         (= c \{) (recur (next stream) (str re \() (inc curly-depth))
         (= c \}) (recur (next stream) (str re \)) (dec curly-depth))
         (and (= c \,) (< 0 curly-depth)) (recur (next stream) (str re \|)
                                                 curly-depth)
         (#{\. \( \) \| \+ \^ \$ \@ \%} c) (recur (next stream) (str re \\ c)
                                                  curly-depth)
         :else (recur (next stream) (str re c) curly-depth)))))

(defn glob [pattern]
  "Returns files matching glob pattern."
  (let [parts (split pattern)
        root (if (= (count parts) 1) "." (apply join (butlast parts)))
        regex (glob->regex (last parts))]
    (map #(.getPath %) (seq (.listFiles (File. root)
                                        (reify FilenameFilter
                                          (accept [_ _ filename]
                                            (if (re-find regex filename)
                                              true false))))))))



; walk helper functions
(defn- w-directory? [f]
  (.isDirectory f))
(defn- w-file? [f]
  (.isFile f))
(defn- w-children [f]
  (.listFiles f))
(defn- w-base [f]
  (.getName f))

; FIXME: I'm sure the Clojure gurus out there will make this a 1 liner :)
(defn walk [path func]
  "Walk over directory structure. Calls 'func' with [root dirs files]"
  (loop [loc (zip/zipper w-directory? w-children nil (io/file path))]
    (when (not (zip/end? loc))
      (let [file (zip/node loc)]
        (if (w-file? file)
          (recur (zip/next loc))
          (let [kids (w-children file)
                dirs (set (map w-base (filter w-directory? kids)))
                files (set (map w-base (filter w-file? kids)))]
            (func (strinfify file) dirs files)
            (recur (zip/next loc))))))))
