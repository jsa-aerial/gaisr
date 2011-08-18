(ns zscore
  (:use edu.bc.utils)
  (:require [clojure.contrib.io :as io]
	    [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]))

(defmacro <=_< [a b c]
  `(let [a# ~a
         b# ~b
         c# ~c]
     (and (<= a# b#) (< b# c#))))

(defn next-base [gcat at gc]
  (let [r (rand)]
    ;;(pr r (* gc gcat) (+ gcat (* at (- 1 gcat))))
    (cond
     (< r (* gc gcat)) "g"
     (<=_< (* gc gcat) r gcat) "c"
     (<=_< gcat r (+ gcat (* at (- 1 gcat)))) "a"
     :else "t")))

(defn generate-seq [length gcat at gc]
  (seq (repeatedly length #(next-base gcat at gc))))


(defn generate-set-file [len filespec]
  (io/with-out-writer filespec
    (doseq [gcat (range 0.25 0.76 0.05)
            at (range 0.25 0.76 0.05)
            gc (range 0.25 0.76 0.05)]
      (println ">" [len gcat at gc])
      (doseq [x (repeatedly 1000 #(str/join "" (generate-seq len gcat at gc)))]
        (println x))))
  filespec)

(defn generate-set-files [dir]
  (pmap #(generate-set-file %1 (str dir "seq-data-" % ".txt"))
        (range 50 401 50)))


(defn ratio [m k1 k2]
  (let [sum (+ (get m k1 0) (get m k2 0))]
    (if (zero? sum) 0 (double (/ (get m k1 0) sum)))))

(defn between? [theory actual]
  (< (- theory 0.01) actual (+ 0.01 theory)))

(defn read-file2 [file]
  (doseq [ln (io/read-lines file)]
    (let [[p s] (str/split #" , " ln)
          param (read-string p)
          seq (drop 1 (str/split #"" s))
          gcat (param 1)
          at (param 2)
          gc (param 3)
          fr-per (reduce (fn [m [base freq]]
                           (assoc m base (double (/ freq (param 0)))))
                         {} (frequencies seq))
          gcat-per (double (+ (get fr-per "g" 0) (get fr-per "c" 0)))
          at-per (ratio fr-per "a" "t")
          gc-per (ratio fr-per "g" "c")
          fr-per (assoc fr-per "gcat" gcat "at" at "gc" gc)]
      (when-not (and (between? at at-per) (between? gc gc-per)
		     (between? gcat gcat-per))
        [p (count s) [at at-per] [gc gc-per] [gcat gcat-per]]))))




(do-text-file ["/home/peis/bin/overnight-set.txt"]
   :x)



(defn ratio [m k1 k2]
  (double (/ (get m k1 0) (+ (get m k1 0) (get m k2 0)))))

(defn stat [in-set]
  (map (fn [x]
         (map (fn [[param y]]
                (let [fr (frequencies (flatten y))
                      total (reduce (fn [s [k v]]
                                      (+ s v))
                                    0 fr)
                      fr-per (reduce (fn [m [k v]]
                                       (assoc m k (double (/ v total))))
                                     {} fr)
                      gcat (double (+ (get fr-per "g" 0) (get fr-per "c" 0)))
                      at (ratio fr-per "a" "t")
                      gc (ratio fr-per "g" "c")
                      fr-per (assoc fr-per "gcat" gcat "at" at "gc" gc)]
                  [param (sort-by first fr-per) total]))
              x))
       in-set))

(defn print-stat [in-set]
  (doseq [i (stat in-set)]
    (doseq [j i]
      (let [[p per n] j
            gcat (nth p 1)
            at (nth p 2)
            gc (nth p 3)
            at-per (second (nth per 1))
            gc-per (second (nth per 4))
            gcat-per (second (nth per 5))
            between? (fn [theory actual] (< (* 0.99 theory) actual (* 1.01 theory))) ]
        (when-not (or (between? at at-per) (between? gc gc-per) (between? gcat gcat-per))
          (prn p ", at" [at at-per] ", gc" [gc gc-per] ", gcat" [gcat gcat-per]))))))


;;(def runovernight
;;  (future (io/with-out-writer "/home/peis/bin/overnight-pmap.txt"
;;           (doseq [i (stat)] (doseq [j i] (prn j))))))

;; (def runovernight2
;;  (future (io/with-out-writer "/home/peis/bin/overnight-set.txt"
;;            (print-set (generate-set)))))
