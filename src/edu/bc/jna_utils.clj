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
            ;;[clojure.contrib.java-utils :as ju]
            ;;[edu.bc.fs :as fs]
            )

  (:import ;;[java.lang.reflect AccessibleObject]
           [com.sun.jna Native Pointer])

  (:use net.n01se.clojure-jna))


;;; The key pieces are the Native and Pointer classes from jna and
;;; clojure-jna which wraps a more dynamic and Clojure like interface
;;; around some of it.
;;;
;;; See http://twall.github.com/jna/3.4.0/javadoc/ for details on JNA


(defn jna-malloc
  "Create a 'C' level USB buffer of size SIZE.  Returns a pair (as a
   vector): [jptr nptr].  Where jptr is a jna Pointer object which
   encapsulates the native pointer nptr and its associated buffer.
   nptr is native C ptr object (a long...).  jptr has all the
   primitive type get and set methods publically available so you can
   manipulate the buffer without having to actually touch it
   directly."
  [size]
   (let [nptr (Native/malloc size)]
     [(Pointer. nptr) nptr]))


(jna-invoke Integer c/printf "Hello there, %s\n" "Clojure")


;;; Test out getting Disk stat information by calling stdlib statvfs.
;;; This requires the C struct statvfs structure - or more to the
;;; point here - memory sufficient to hold it.  This is not obvious
;;; (and at first this was 88 and was blowing things up).  It's
;;; recommended on things like this to get the true sizes from a quick
;;; C program (via sizeof) and then use them.  Not clear/sure about a
;;; more automated variation of that which could be all done at
;;; runtime.
;;;
(let [[jptr nptr] (jna-malloc 112)]
  (println nptr)
  (jna-invoke Integer c/statvfs "/data2" nptr)
  (println nptr)
  ;; Be wary/careful of machine word size!!!
  (let [fbsize (. jptr getLong 0)
        frsize (. jptr getLong 8)
        blocks (. jptr getInt 16)
        bfree  (. jptr getInt 24)
        bavail (. jptr getInt 32)
        files  (. jptr getInt 40)
        ffree  (. jptr getInt 48)
        favail (. jptr getInt 56)]

    (println nptr)

    (println "f_fbsize" fbsize)
    (println "f_frsize" frsize)
    (println "blocks" (* 4 blocks))
    (println "bfree" (* 4 bfree))
    (println "bavail" (* 4 bavail))
    (println "files" (* 4 files))
    (println "ffree" (* 4 ffree))
    (println "favail" (* 4 favail))

    (Native/free nptr)
    ))



;;; Test making sure that free, frees our malloc'd memory.  That is,
;;; the following does not leak.  NOTE: jna-malloc'd memory must be
;;; Native/free'd as it is outside the Java heap and the GC doesn't
;;; know about it!
;;;
(defn test-free [n]
  (let [stg "CGCAGGGAUACCCGCG"]
    (dotimes [i n]
      (let [[jptr nptr] (jna-malloc (inc (count stg)))
            e (jna-invoke Float RNA/fold stg nptr)]
        (Native/free nptr)
        ))))




;;; Let's call RNAfold directly!  Must build a sharable lib from the
;;; Vienna .os first and place in library path.  Here it is libRNA.so
;;; Watch out for version 2.0 - it is inconsistent (undefined symbols
;;; referenced!), but you can hack it to work.
;;;
(let [stg "CGCAGGGAUACCCGCG"
      [jptr nptr] (jna-malloc (inc (count stg)))
      e (jna-invoke Float RNA/fold stg nptr)]
  (println jptr :===> nptr)
  (println "E =" e)
  (println "Structure =" (. jptr getString 0 false))
  (Native/free nptr))

(jna-invoke Float RNA/energy_of_struct "CGCAGGGAUACCCGCG" "(((.((....)).)))")


;;; Alifold and energy_of_alistruct....
;;;
(let [stg "CGCAGGGAUACCCGCG"
      [in-jptr in-nptr] (jna-malloc (inc (count stg)))
      [out-jptr out-nptr] (jna-malloc (inc (count stg)))]
  (. in-jptr setString 0 stg false)
  (println "B4 Call" (. in-jptr getString 0 false))
  (println (jna-invoke Float RNA/fold in-nptr out-nptr))
  (println "Structure =" (. out-jptr getString 0 false))
  (Native/free in-nptr)
  (Native/free out-nptr)
  )


(let [stgarray (into-array
                ["AGAUAUAUAUGCUAAGCGCCGCAGACAGCGGGU.GC....GUU....U.G"
                 "CCCCAGUCGAAUAGCAAGCCGAAGACAGCAGGU.GCC..CGCG.....GG"
                 "AAUAAAAUAACAUUUGUACCGUAGACAGCAGGU.GC....GAC....A.G"
                 "GUUUUGAUAAAAAACUGACCUAAGACAGCAGGG.GAG..CAUG.....CU"])
      ;; structure space to be filled in with structure
      [out-jptr out-nptr] (jna-malloc (inc (count (first stgarray))))
      ;; Need array of floats for output energies
      [eout-jptr eout-nptr] (jna-malloc (* 2 4))] ; 2 * 4 bytes per float
  ;; Amazingly, JNA just knows how to correctly pass an array of
  ;; strings as a char**
  (println (jna-invoke Float RNA/alifold stgarray out-nptr))
  (println "Structure =" (. out-jptr getString 0 false))
  ;;
  (println "Alistruct free engery ret "
           (jna-invoke Float RNA/energy_of_alistruct
                       stgarray out-nptr (count stgarray) eout-nptr))
  (println "energies" (. eout-jptr getFloat 0) (. eout-jptr getFloat 4))
  (println "Structure =" (. out-jptr getString 0 false)))




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
(let [[fa-jptr fa-nptr] (jna-malloc (* 8 6))
      [ia-jptr ia-nptr] (jna-malloc (* 4 2))
      fa (double-array [1.0 2.0 3.0 4.0 5.0])]

  ;; Put the Doubles into the "array"
  (doseq [i (take 5 (iterate inc 0))]
    (let [f (aget fa i)]
      (println "f =" f "i =" i)
      (. fa-jptr setDouble (* 8 i) f)))

  ;; Set the integer count of array "size"
  (. ia-jptr setInt 0 5)

  ;; Are they in the "arrays" OK??
  (doseq [i [0 1 2 3 4]]
    (println (. fa-jptr getDouble (* 8 i))))
  (println (. ia-jptr getInt 0))

  ;; Call foo and print return
  (println (jna-invoke Void RTest/foo ia-nptr fa-nptr))

  ;; Check to see if results are there.
  (doseq [i [0 1 2 3 4]]
    (println (. fa-jptr getDouble (* 8 i))))

  (Native/free fa-nptr)
  (Native/free ia-nptr))

