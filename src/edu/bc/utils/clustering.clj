;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                  U T I L S . C L U S T E R I N G                         ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011-2013 Trustees of Boston College                       ;;
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

(ns edu.bc.utils.clustering

  "Various data clustering algorithms, techniques, functions.
   Generally applicable, but typically used for sequence clustering of
   various sorts."

  (:require [clojure.contrib.math :as math]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        edu.bc.utils.probs-stats
        [clojure.pprint
         :only [cl-format]]))




(defn edist
  "Simple Euclidean distance between points x and y.  x and y are
   points from the same dimensional space (1-d, 2-d, 3-d, ... n-d)
  "
  [x y]
  (cond
   (and (number? x) (number? y))
   (math/abs (- x y))
   (and (vector? x) (vector? y))
   (vecdist x y)
   :else
   (raise :type :illegal-params
          :msg "edist must have two numbers or vectors")))


(defn extreme-pd
  "For point X, compute its 'extreme' point and its corresponding
   distance from all the points in COLL.  The 'extremeness' is
   determined by COMP a comparison predicate.  For example, comp = <
   indicates 'nearest'.

   DISTFN is the distance function used to compute the distances of x
   from any p in coll.

   If coll is a map, uses (vals coll) for the candidate collection of
   points.
  "
  [x distfn comp coll]
  (let [coll (if (map? coll) (vals coll) coll)
        pds (map #(do [% (distfn % x)]) coll)]
    (reduce (fn[[mc md :as m] [c d :as p]] (if (comp d md) p m))
            (first pds) (rest pds))))

(defn nearest-pd
  "Return the nearest point p in coll to x along with its distance.
   Distance is given by distfn.  If coll is a map uses (vals coll).
   Returns [p,d], where p is the nearest point and d its distance.
  "
  [x distfn coll]
  (extreme-pd x distfn < coll))

(defn nearest
  "Return the nearest point p in coll to x.  Distance is given by
   distfn.  If coll is a map uses (vals coll).
  "
  [x distfn coll]
  (first (nearest-pd x distfn coll)))

(defn farthest-pd
  "Return the farthest point p in coll to x along with its distance.
   Distance is given by distfn.  If coll is a map uses (vals coll).
   Returns [p,d], where p is the farthest point and d its distance.
  "
  [x distfn coll]
  (extreme-pd x distfn > coll))

(defn farthest
  "Return the farthest point p in coll to x.  Distance is given by
   distfn.  If coll is a map uses (vals coll).
  "
  [x distfn coll]
  (first (farthest-pd x distfn coll)))


(defn ith-sum-sqr-err
  "Sum of Squared Error for cluster Ci.  Ci is a map entry with key
   the mean mi of Ci and val the points assigned to Ci.  distfn is the
   distance function for the err (and must be the same as the distance
   function assigning the points to Ci).
  "
  [distfn Ci]
  (let [[mi xs] Ci]
    (sum #(sqr (distfn mi %)) xs)))

(defn sum-sqr-err
  "Sum of Squared Error for a set of clusters (result of single loyd
   step or end of k-means or ...).  distfn is the distance function
   used to form the clusters (in a loyd step) and is also the distance
   function for computing the errors.
   "
  [distfn clusters]
  (sum (partial ith-sum-sqr-err distfn) clusters))

(defn centers
  "Compute a new set of centers from CLUSTERS using AVGFN as a 'mean'
   for the points in each cluster Ci of clusters.  Returns an 'eager'
   seq of these new centers.
  "
  [avgfn clusters]
  (reduce (fn[S cl] (conj S (avgfn cl)))
          [] clusters))

(defn clusters
  "Form and return a set of clusters.  Each cluster is the set of
   points from COLL closest to a point ci in CENTERS as determined by
   distance function distfn (see nearest).  Result is a map with keys
   ci and values the cluster points for ci.
  "
  [distfn coll centers]
  (reduce (fn[M [x c]] (assoc M c (conj (get M c []) x)))
          {} (map #(do [% (nearest % distfn centers)]) coll)))


(defn split-worst-cluster
  "Clusters a cluster map as given by function CLUSTER (for which
   see).  The functions distfn and avgfn are the distance and 'mean'
   functions used to build clusters.

   Clusters is assumed to be smaller in count than it 'should' be.
   Using an error measure, which here is SSE (see ith-sum-sqr-err and
   sum-sqr-err), clusters has a worst cluster Ci.  Split Ci by first
   sorting its points by their distance from cmi, the center of Ci.
   Then remove the points within WRADIUS percent of the furthest point
   and place in new cluster wCi.  Remove Ci from clusters and add the
   pair of clusters from splitting Ci:

   new-clusters = clusters - Ci + (Ci - wCi) + wCi

   Return new-clusters.
  "
  [distfn avgfn clusters & {:keys [wradius] :or {wradius 0.8}}]
  (let [vclus (ensure-vec clusters)
        errs (sort-by second > (map #(do [%1 (ith-sum-sqr-err distfn %2)])
                                    (iterate inc 0) vclus))
        worst2best (map #(vclus(first %)) errs)
        [wmi wxs] (first (drop-until #(> (count (second %)) 1) worst2best))
        newclusters (dissoc clusters wmi)
        ranked-wxs (sort-by #(distfn wmi %) > wxs)
        wscore (* wradius (distfn wmi (first ranked-wxs)))
        [nc1 nc2] [(take-while #(> (distfn wmi %) wscore) ranked-wxs)
                   (drop-until #(<= (distfn wmi %) wscore) ranked-wxs)]]
    (when (not (seq ranked-wxs))
      (raise :type :unresolvable-cluster
             :fn 'split-worst-cluster
             :args [distfn clusters]
             :errs errs
             :worst2best worst2best))
    (assoc newclusters (avgfn nc1) nc1 (avgfn nc2) nc2)))


(defn loyd-step
  "Primary step in kmeans algorithm.  Recluster data to the 'new'
   centers incenters.  Take new clusters and compute and return newer
   centers for all new clusters.  That's the basic Loyd step, but here
   there is an extra wrinkle:

   If (count new-clusters) is not equal (count incenters), some input
   clusters have coalesced and need to be resplit (see
   split-worst-cluster).  The result is then 'restepped' until the
   count of new clusters is count of incenters.
  "
  [distfn avgfn data incenters]
  (loop [clus (clusters distfn data incenters)]
    (if (= (count clus) (count incenters))
      (->> clus vals (centers avgfn))
      (recur (split-worst-cluster distfn avgfn clus)))))


(defn kmeans
  "Computes kmeans clustering over COLL starting with initial set of
   centers INITIAL-CENTERS, using DISTFN as the distance function
   between points and AVGFN as the 'means' function for a cluster.

   let k = count initial-centers
   repeat
     form k clusters by grouping all p in coll with nearest center
     compute new k centers from clusters
   until centers-i = centers-i+1

   Return [Cs sse], where Cs is the set of clusters asssociated with
   final centers and sse is the sum-sqr-err of Cs
  "
  [initial-centers coll
   & {:keys [distfn avgfn] :or {distfn vecdist avgfn vecmean}}]
  (->> (iterate (partial loyd-step distfn avgfn coll) initial-centers)
       take-until-nochange
       last
       (clusters distfn coll)
       (#(do [% (sum-sqr-err distfn %)]))))


;;; (km++init 3 [0, 1, 2, 3, 4] edist)
(defn- km++init
  "The kmeans++ initialization computation for K clusters over data
   COLL using distance measure DISTFN:

   let D2i (fn[] (apply min (map (sqr (distfn ci %)) C)))
       C = {c1}, c1 random uniformly picked from coll
   loop [C C]
     if (= (count C) k)
       C
       let ci (pick ci from coll with prob (/ (D2i ci) (sum (D2i %) C)))
         recur (conj C ci)

   This weights the probability of picking ci by its distance from
   current centers and so, those further away from all current centers
   have higer chance of being picked.  However, outliers have lower
   overall chance than something from _groups_ positioned further away
   as the _sum_ of their probabilities swamps the outlier.  So, tends
   to maximize initial centers from all natural k clusters, thereby
   ensuring:

   1. much faster convergence
   2. convergence to clusters much closer to real clusters
  "
  [k distfn coll]
  (let [coll (ensure-vec coll)]
    (loop [Cs [(rand-nth coll)]]
      (if (= (count Cs) k)
        Cs
        (let [D2s (map (fn[x] (apply min (map #(sqr (distfn % x)) Cs))) coll)
              Dtot (sum D2s)
              prs (map #(/ % Dtot) D2s)
              r (rand)
              i (ffirst (drop-while
                         (fn[[i p]] (> r p))
                         (map #(do [%1 %2])
                              (iterate inc 0)
                              (cumsum prs))))]
          (recur (conj Cs (coll i))))))))

(defn kmeans++
  "kmeans++ clustering of COLL into K clusters.  Picks initial-centers
   as (k++init k distfn coll) and then proceeds via kmeans.  As for
   kmeans, DISTFN is the distance function between points and AVGFN is
   the 'means' function for a cluster.
  "
  [k coll & {:keys [distfn avgfn] :or {distfn vecdist avgfn vecmean}}]
  (kmeans (km++init k distfn coll) coll :distfn distfn :avgfn avgfn))





;;; ----------------- Ad Hoc Testing Stuff -------------------------------

(comment

(def data1 [2 3 5 6 10 11 100 101 102])
(def data (for [i (range 50)] (int (rand 100))))

(->> (iterate (partial loyd-step edist #(double (/ (sum %) (count %))) data)
              (km++init 6 edist data))
     take-until-nochange
     (map #(clusters edist data %))
     (map #(sum-sqr-err edist %))
     sort)

(clusters edist data1 [0 10])
(clusters edist (range 5) (loyd-step edist mean (range 5) [0 1]))

(->> (iterate (partial loyd-step edist mean data1) [1 2 3])
     take-until-nochange
     (map #(clusters edist data1 %))
     last
     (sum-sqr-err edist))

(let [data (range 5)]
  (map #(vals (clusters edist data %))
       (take-until-nochange
        (iterate (partial loyd-step edist mean data)
                 [0 1]))))

(let [data [[1 2 3] [3 2 1] [100 200 300] [300 200 100] [50 50 50]]]
  (->> (iterate (partial loyd-step vecdist vecmean data)
                [[1 1 1] [2 2 2] [3 3 3]])
       take-until-nochange
       (map #(clusters vecdist data %))))

(let [data [[1 2 3] [3 2 1] [100 200 300] [300 200 100] [50 50 50]]]
  (for [k (range (dec (count data)) 0 -1)]
    [k (kmeans++ k data)]))


)