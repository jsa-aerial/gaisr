(ns edu.bc.jna-utils

  "Interfacing with native C lib stuff.  This includes all manner of
   things from stuff like stdlib, statvfs, string lib, to rnafold native,
   libsvm, etc.
   *** CURRENTLY JUST PROOF OF CONCEPT AND EXAMPLES ***"

  (:require [clojure.contrib.math :as math]
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.contrib.io :as io]
            [clojure.contrib.properties :as prop]
            [edu.bc.fs :as fs])

  (:use [clojure.contrib.condition
         :only [raise handler-case *condition*
                print-stack-trace stack-trace-info]]

        [clojure.contrib.pprint
         :only (cl-format compile-format)]

        net.n01se.clojure-jna)) ; The key piece


(defn jna-malloc
  "Create a 'C' level USB buffer of size SIZE.  Returns a pair (as a
   vector): [ptr-to-the-buffer the-buffer] Where ptr-to-the-buffer is
   a JNA/C ptr object and the-buffer is a java.nio.DirectByteBuffer
   object." [size]
  (let [buf (make-cbuf size)
        ptr (pointer buf)]
    [ptr buf]))


(jna-invoke Integer c/printf "Hello there, %s\n" "Clojure")

(let [[ptr buf] (jna-malloc 88)]
  (jna-invoke Integer c/statvfs "/data2" ptr)
  ;; Be wary/careful of machine word size!!!
  (let [fbsize (.getLong buf)
        frsize (.getLong buf 8)
        blocks (.getInt buf 16)
        bfree  (.getInt buf 24)
        bavail (.getInt buf 32)
        files  (.getInt buf 40)
        ffree  (.getInt buf 48)
        favail (.getInt buf 56)]


    (println "f_fbsize" fbsize)
    (println "f_frsize" frsize)
    (println "blocks" (* 4 blocks))
    (println "bfree" (* 4 bfree))
    (println "bavail" (* 4 bavail))
    (println "files" (* 4 files))
    (println "ffree" (* 4 ffree))
    (println "favail" (* 4 favail))
    ))

;;; Let's call RNAfold directly!  Must build a sharable lib from the
;;; Vienna .os first and place in library path.  Here it is libRNA.so
;;; Watch out for version 2.0 - it is inconsistent (undefined symbols
;;; referenced!), but you can hack it to work.
;;;
(let [stg "CGCAGGGAUACCCGCG"
      [ptr buf] (jna-malloc (inc (count stg)))
      e (jna-invoke Float RNA/fold stg ptr)]
  (println "E =" e)
  (println "Structure =" (.getString ptr 0 false)))

(jna-invoke Float RNA/energy_of_struct "CGCAGGGAUACCCGCG" "(((.((....)).)))")


;;; Alifold and energy_of_alistruct....
;;;
(let [stg "CGCAGGGAUACCCGCG"
      [inptr inbuf] (jna-malloc (inc (count stg)))
      [outptr outbuf] (jna-malloc (inc (count stg)))]
  (.setString inptr 0 stg false)
  (println "B4 Call" (.getString inptr 0 false))
  (println (jna-invoke Float RNA/fold inptr outptr))
  (println "Structure =" (.getString outptr 0 false))
  )

(let [stgarray (into-array
                ["AGAUAUAUAUGCUAAGCGCCGCAGACAGCGGGU.GC....GUU....U.G"
                 "CCCCAGUCGAAUAGCAAGCCGAAGACAGCAGGU.GCC..CGCG.....GG"
                 "AAUAAAAUAACAUUUGUACCGUAGACAGCAGGU.GC....GAC....A.G"
                 "GUUUUGAUAAAAAACUGACCUAAGACAGCAGGG.GAG..CAUG.....CU"])
      ;; structure space to be filled in with structure
      [outptr outbuf] (jna-malloc (inc (count (first stgarray))))
      ;; Need array of floats for output energies
      [eoutptr eoutbuf] (jna-malloc (* 2 4))] ; 2 * 4 bytes per float
  ;; Amazingly, JNA just knows how to correctly pass an array of
  ;; strings as a char**
  (println (jna-invoke Float RNA/alifold stgarray outptr))
  (println "Structure =" (.getString outptr 0 false))
  ;;
  (println "Alistruct free engery ret "
           (jna-invoke Float RNA/energy_of_alistruct
                       stgarray outptr (count stgarray) eoutptr))
  (println "energies" (.getFloat eoutbuf) (.getFloat eoutbuf 4))
  (println "Structure =" (.getString outptr 0 false)))




;;; Example 'arbitrary' C func.
;;;
;;; void foo(int *nin, double *x)
;;; {
;;;     int n = nin[0];
;;;
;;;     int i;
;;;
;;;     for (i=0; i<n; i++)
;;;             x[i] = x[i] * x[i];
;;; }
;;;
;;; gcc -c -fPIC foo.c -o foo.o
;;; gcc -shared -Wl,-soname,libRTest.so -o libRTest.so  foo.o
;;;
;;; Everything seems to work off C allocated "array" things (makes
;;; sense since that is what C does)
;;;
(let [[faptr fabuf] (jna-malloc (* 8 6))
      [iaptr iabuf] (jna-malloc (* 4 2))
      fa (double-array [1.0 2.0 3.0 4.0 5.0])]

  ;; Put the Doubles into the "array"
  (doseq [i (take 5 (iterate inc 0))]
    (let [f (aget fa i)]
      (println "f =" f "i =" i)
      (.setDouble faptr (* 8 i) f)))

  ;; Set the integer count of array "size"
  (.setInt iaptr 0 5)

  ;; Are they in the "arrays" OK??
  (doseq [i [0 1 2 3 4]]
    (println (.getDouble fabuf (* 8 i))))
  (println (.getInt iabuf))

  ;; Call foo and print return
  (println (jna-invoke Void RTest/foo iaptr faptr))

  ;; Check to see if results are there.
  (doseq [i [0 1 2 3 4]]
    (println (.getDouble fabuf (* 8 i)))))
