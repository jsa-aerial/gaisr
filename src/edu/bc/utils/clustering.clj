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
            [clojure.set :as set]
            [edu.bc.fs :as fs]
            [edu.bc.utils.graphs :as gr])

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


(defn dist-matrix
  "Computes and returns the pairwise distance matrix for items in
   COLL.  The distances between items is given by distfn.  SYM
   indicates that distfn is symmetric (default) and keyfn returns a
   key suitable for map entries for an item and defaults to identity.
   The idea behind keyfn is that some items may be large or otherwise
   complicated and thus expensive to compare but a unique key may be
   generated for them beforehand.  A typical scenario would be to
   'Goedel number' them.

   Returns {[k v] | k=[(kefn i) (kefn j)] v =(distfn i j)}

   Uses reducers and xfold to parallelize computation with auto
   computed queue granularity (see xfold)
  "
   [distfn coll & {:keys [sym keyfn] :or {sym true keyfn identity}}]
   (let [coll (vec coll)
         cnt (count coll)]
     (reduce (fn[M [ik jk d]]
               (if sym
                 (assoc M [ik jk] d [jk ik] d)
                 (assoc M [ik jk] d)))
             {} (xfold (fn [[i j]] [(keyfn i) (keyfn j) (distfn i j)])
                       (for [k (range cnt)
                             l (range cnt)
                             :let [i (coll k)
                                   j (coll l)]
                             :when (if sym (< k l) (not= k l))]
                         [i j])))))


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
  "Sum of Squares Error for cluster Ci.  Ci is a map entry with key
   the mean mi of Ci and val the points assigned to Ci.  distfn is the
   distance function for the err (and must be the same as the distance
   function assigning the points to Ci).
  "
  [distfn Ci]
  (let [[mi xs] Ci]
    (sum #(sqr (distfn mi %)) xs)))

(defn sum-sqr-err
  "Sum of Squares Error for a set of clusters (result of single loyd
   step or end of k-means or ...).  distfn is the distance function
   used to form the clusters (in a loyd step) and is also the distance
   function for computing the errors.
   "
  [distfn clusters]
  (sum (partial ith-sum-sqr-err distfn) clusters))



;;; ------------------------------------------------------------------------;;;
;;;                   Loyd step and flavors of kmeans                       ;;;

(defn centers
  "Compute a new set of centers from CLUSTERS using AVGFN as a 'mean'
   for the points in each cluster Ci of clusters.  Returns an 'eager'
   seq of these new centers.
  "
  [avgfn clusters]
  (xfold (fn[cl] (avgfn cl)) clusters))

(defn clusters
  "Form and return a set of clusters.  Each cluster is the set of
   points from COLL closest to a point ci in CENTERS as determined by
   distance function distfn (see nearest).  Result is a map with keys
   ci and values the cluster points for ci.
  "
  [distfn coll centers]
  (reduce (fn[M [x c]] (assoc M c (conj (get M c []) x)))
          {} (xfold #(do [% (nearest % distfn centers)]) coll)))


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
        worst2best (map #(vclus (first %)) errs)
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
       (clusters distfn coll)))


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
  (let [coll (ensure-vec (set coll))]
    (loop [Cs #{(rand-nth coll)}]
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



;;; ------------------------------------------------------------------------;;;
;;;         KNN and KRNN based clustering (RECORD, CHAMELEON)               ;;;


(defn knn
  "Compute K nearest neighbors of point P in collection COLL.
   Distance is given by means of distfn, preferably a metric, but can
   be a 'similarity measure' such as relative entropy.  O(nlogn)
   complexity.
  "
  [k distfn p
  coll]
  (take k (sort-by #(distfn p %) coll)))

(defn knn-graph
  "Compute the k nearest neighbors graph over the items in COLL as
   determined by distance function distfn, preferably a metric, but
   can be a 'similarity measure' such as relative entropy.  Returns
   the graph encoded as a map (sparse edge set matrix), with keys
   items in COLL and values k-sets of items in COLL.
  "
  [k distfn coll]
  (let [coll (set coll)]
    (reduce
     (fn[G p] (assoc G p (knn k distfn p (set/difference coll #{p}))))
     {} coll)))


(defn krnn-graph
  "Compute the k reverse nearest neighbors graph over the items in
   COLL as determined by the k nearest neighbors graph over COLL (see
   knn-graph).  Distance of knn graph is determined by distance
   function distfn, preferably a metric, but can be a 'similarity
   measure' such as relative entropy.  Returns the graph encoded as a
   map (sparse edge set matrix), with keys items in COLL and values
   sets of items in COLL.

   NOTE: the size of the value sets for krnn need not be k, typically
   isn't k, and can range from 0 to (count coll).
  "
  [k distfn coll]
  (let [knngrph (knn-graph k distfn coll)]
    (conj (reduce (fn[[krnnM rnncntM] p]
                    (reduce (fn[[krnnM rnncntM] q]
                              [(assoc krnnM q (conj (get krnnM q []) p))
                               (assoc rnncntM q (inc (get rnncntM q 0)))])
                            [krnnM rnncntM] (knngrph p)))
                  [{} {}] coll)
          knngrph)))


(defn split-krnn
  "Takes a krnn graph (as built by krnn-graph) for value k, with point
   counts given by rnncntM (as built by krnn-graph) and the
   corresponding knngrph that is the basis of the krnn graph, and
   returns two new graphs [rnnG>k rnnG<k]:

   rnnG>k is the subgraph of krnngrph whose nodes have >= k edges

   rnnG<k is the subgraph of krnngrph whose nodes have < k edges

   Additionally, all edge set nodes in both subgraphs that do not
   appear as nodes in the subgraphs are removed (another choice would
   have been to include them with null edge sets, but that would have
   violated the constraints on rnnG>k).

   These two graphs form the basis for initial sets of clusters based
   on the set of strongly connected components (SCC) in them.  These
   SCC form the basis of both the Chameleon and RECORD clustering
   algorithms.
  "
  [k krnngrph rnncntM knngrph]
  (let [[rnnG>k rnnG<k]
        (reduce (fn[[rnnG>k rnnG<k] p]
                  (if (< (rnncntM p) k)
                    (let [rnnG>k (reduce
                                  (fn[G q]
                                    (if (G q)
                                      (assoc G q (remove #(= % p) (G q)))
                                      G)) ; If it was removed, don't put back!
                                  rnnG>k (knngrph p))]
                      [(dissoc rnnG>k p)
                       (assoc rnnG<k p (rnnG>k p))])
                    [rnnG>k rnnG<k]))
                [krnngrph {}] (keys rnncntM))
        ;; Make final pass on just constructed outliers removing edge
        ;; set points in rnnG>k (well mostly - really any not in
        ;; rnnG<k)
        rnnG<k (let [nodes (keys rnnG<k)]
                 (reduce (fn[rnnG<k q]
                           (if (not (in q nodes))
                             (reduce (fn[G n]
                                       (assoc G n (remove #(= % q) (G n))))
                                     rnnG<k nodes)
                             rnnG<k))
                         rnnG<k (-> rnnG<k vals flatten set)))]
    [rnnG>k rnnG<k]))


(defn refoldin-outliers
  "Takes a krnn graph and the SCC 'clusters' corresponding to the s"
  [krnngrph rnnG>k-sccs rnnG<k-sccs]
  (let [clusters (first (reduce (fn[[M i] scc]
                                  [(assoc M i scc) (inc i)])
                                [{} 0] rnnG>k-sccs))
        outliers (->> rnnG<k-sccs
                      (map seq) flatten
                      (map #(do [% (set (krnngrph %))]))
                      (into {}))]
    #_(prn clusters)
    (->>
     (reduce (fn[clus [o ons]]
               (->> clus (map (fn[[k v]] [k (count (set/intersection v ons))]))
                    (sort-by second >) first
                    ((fn[[k cnt]]
                       (if (not= 0 cnt)
                         (assoc clus k (conj (clus k) o))
                         (assoc clus (gen-uid) #{o}))))))
             clusters outliers)
     vals set)))


(defn krnn-clust
  ""
  [k distfn coll & {:keys [keyfn] :or {keyfn identity}}]
  (let [dm (dist-matrix distfn coll :keyfn keyfn)
        kcoll (map keyfn coll)
        [krnngrph rnncntM knngrph] (krnn-graph k #(get dm [%1 %2]) kcoll)]
    (->> (split-krnn k krnngrph rnncntM knngrph)
         (map #(gr/tarjan (keys %) %))
         (apply refoldin-outliers krnngrph))))




;;; --------------- Clustering Distance and Validity Measures --------------;;;


(defn center-distances
  "Pairwise distances of centers of clusters Ci in CLUSTERING by means
   of distance function DISTFN.
  "
  [distfn clustering]
  (map (fn[[[ci _] [cj _]]] (distfn ci cj)) (combins 2 clustering)))

(defn center-dist-expect
  "The expectation of all center distances of clusters in CLUSTERING.
   Distances are computed by DISTFN.  Expectation is computed by AVGFN
   which defaults to the mean.
  "
  [distfn clustering & {:keys [avgfn] :or {avgfn mean}}]
  (avgfn (center-distances distfn clustering)))


(defn intra-distances
  "Pairwise distances of points in a cluster CLUSTER as given by
   DISTFN.  CLUSTER is either an element of a result of (clusters
   ...), i.e., a map entry with center key and points val, or the
   point collection of such cluster.
  "
  [distfn cluster]
  (let [clu (if (map-entry? cluster) (val cluster) cluster)]
    (map (fn[[x y]] (distfn x y)) (combins 2 clu))))

(defn intra-dist-expect
  "The expectation of all pairwise point distances in cluster CLUSTER.
   Distances are computed by DISTFN and expectation by AVGFN, which
   defaults to mean. CLUSTER is either an element of a result
   of (clusters ...), i.e., a map entry with center key and points
   val, or the point collection of such cluster.
  "
  [distfn cluster & {:keys [avgfn] :or {avgfn mean}}]
  (avgfn (intra-distances distfn cluster)))


(defn cluster-distances
  "Given a clustering CLUSTERS and its distance function DISTFN,
   returns a set of scores/measures of the quality of the clustering,
   a map with keys [:global, :each, :inter-m-xs, inter-xs, ...].
   Where,

   * global gives a global intra cluster cohesion, uses SSE

   * each gives intra cluster cohesion for each cluster, uses SSE

   * inter a global inter cluster cohesion (sum of sqrs of center
     distances)

   * inter-ms gives pairwise inter cluster cohesion (center distances)

   * inter-m-xs gives pairwise inter cluster cohesion by min distances
     over points xis in Ci to center mj of Cj and vice versa.

   * inter-xs gives pairwise inter cluster cohesion by min distance
     over all [xi, xj] pairs, xi in Ci, xj in Cj
  "
  [distfn clusters]
  (let [cpairs (vec (combins 2 clusters))
        ms (sort-by second
                    (map (fn[l [[mi xis] [mj xjs]]] [l (distfn mi mj)])
                         (iterate inc 0) cpairs))
        mDs (sort-by second
                     (map (fn[l [[mi xis] [mj xjs]]]
                            [l
                             (apply min (map #(distfn mi %) xjs))
                             (apply min (map #(distfn mj %) xis))])
                          (iterate inc 0) cpairs))
        xDs (sort-by second
                     (map (fn[l [[mi xis] [mj xjs]]]
                            [l (apply min (reducem #(do [(distfn %1 %2)])
                                                   concat
                                                   xis xjs))])
                          (iterate inc 0) cpairs))]
    {:global (sum-sqr-err distfn clusters)
     :each (map #(do [((partial ith-sum-sqr-err distfn) %) %]) clusters)
     :inter (sum #(sqr (second %)) ms)
     :inter-ms (map (fn[[i d]] [d (cpairs i)]) ms)
     :inter-m-xs (map (fn[[i d1 d2]] [d1 d2 (cpairs i)]) mDs)
     :inter-xs (map (fn[[i d]] [d (cpairs i)]) xDs)}))


(defn DBI-Si
  "The Davies Bouldin Index of the compactness (aka 'scatter' or Si)
   of cluster CLUSTER.

   Si is the expectation of pairwise distances of points in CLUSTER
   with its center ci.  Distances are given by DISTFN. The expectation
   is AVGFN of these distances.  AVGFN defaults to mean:

   Si = (avgnf (sum (fn[xj] (distfn ci xj)) (val cluster)))

   CLUSTER must be an element of a result of (clusters ...), i.e., a
   map entry with center key and points val.
  "
  [distfn cluster & {:keys [avgfn] :or {avgfn mean}}]
  (let [ci (key cluster)]
    ;; Penalize singleton clusters.  Bit of a hack, would be better
    ;; overall if took the largest so far, but that requires global
    ;; knowledge.  Maybe DBI-Rij could do this???
    (if (< (count (val cluster)) 2)
      (double 10.0)
      #_(math/sqrt (avgfn (map (fn[xj] (sqr (distfn xj ci))) (val cluster))))
      (avgfn (map (fn[xj] (distfn xj ci)) (val cluster)))
      )))

(defn DBI-Rij
  "Rij is a similarity measure between two clusters of a clustering
   based on DBI-Si compactness and center point distance dij.  The
   important aspects are that it is postive definite and symmetric.
   Additionally, 'triangle like' aspects hold:

   * if Sj > Sk and dij = dik, then Rij > Rik

   * if Sj = Sk and dij < dik, then Rij > Rik

   Rij = (/ (+ Si Sj) (distfn ci cj))

   Returns the collection of Rij for all pairwise clusters of
   clustering as a map indexed by i in (range n), n = (count
   clusterings).  The value of an entry is the n-1 set of Ri,i/=j for
   cluster Ci.  Note since Rij=Rji there are duplicate values across
   the entries, but they are computed only once.
  "
  [distfn clustering  & {:keys [avgfn] :or {avgfn mean}}]
  (let [indices (combins 2 (range (count clustering)))
        vclus (vec clustering)]
    (reduce (fn[M [i j :as p]]
              (let [[Ci Cj] [(vclus i) (vclus j)]
                    [ci cj] [(first Ci) (first Cj)]
                    Rij (/ (+ (DBI-Si distfn Ci) (DBI-Si distfn Cj))
                           (distfn ci cj))
                    M (assoc M i (conj (get M i []) Rij))
                    M (assoc M j (conj (get M j []) Rij))]
                M))
            {} indices)))

(defn davies-bouldin-index
  "Computes the Davies Bouldin Index of cluster validity.  This is a
   cluster validity measure which is the average of the maximal Rij
   over all not equal pairwise clusters (see DBI-Rij).
  "
  [distfn clustering & {:keys
  [avgfn] :or {avgfn mean}}]
  (let [Ds (map (fn[[i Rijs]] (apply max Rijs)) ;(median Rijs))
                (DBI-Rij distfn clustering :avgfn avgfn))]
    (mean Ds)))


(defn- ||x|| [v]
  (if (not (coll? v)) v (norm v)))

(defn cluster-stdev
  "Compute the average standard deviation of the spread of clusters in
   CLUSTERING, the result of a (clusters ...) call.  This is not as
   obvious as it may see, and is basically (avg-std-deviation
   clustering) except we cheat a bit as we already have the
   means (centers).
  "
  [distfn avgfn clustering]
  (math/sqrt
   (/ (sum
       (xfold
        (fn[[mi xis]]
          (* (dec (count xis))
             (variance xis :distfn distfn :avgfn avgfn :m mi)))
        1 clustering))
      (- (sum count clustering) (count clustering)))))

(defn density
  "Computes the 'density' of points in COLL relative to u (typically a
   'center' or 'mid point' of some sort) as determined by counts
   within 1 STDEV of u.  Distances are given by distance function
   DISTFN.
  "
  [distfn stdev u coll]
  (let [f (fn[x mi] (if (<= (distfn x mi) stdev) 1 0))]
    (sum (xfold #(f % u) coll))))

(defn intercluster-density
  "Compute the inter cluster point density.  This is the density
   between clusters as determined by point counts 1 cluster-stdev from
   midpoints between cluster centers.

   Note: if (= 1 (count clustering)), simply returns the density of
   the single cluster.

   Returns floating point density score of separation - the smaller
   the better.
  "
  [distfn avgfn clustering & {:keys [stdev]}]
  (let [stdev (if stdev stdev (cluster-stdev distfn avgfn clustering))
        c (count clustering)]
    #_(prn :stdev stdev)
    (if (= c 1)
      (let [[m1 x1s] (first clustering)]
        (density distfn stdev m1 x1s))
      (/ (sum (fn[[mi xis]]
                (sum (fn[[mj xjs]]
                       (if (= mi mj)
                         0
                         (let [u (avgfn mi mj)
                               C (set/union (set xis) (set xjs))]
                           (/ (density distfn stdev u C)
                              (max (density distfn stdev mi xis)
                                   (density distfn stdev mj xjs)
                                   1)))))
                     clustering))
              clustering)
         (* c (dec c))))))

(defn scatt
  "Compute the compactness of a clustering by means of 'density'
   measure.  This is the averaged ratio of variance of clusters in
   clustering to the overall variance of the data set:

   let S all data points
       n (count clustering)
     (* 1/n (sum (fn[Ci] (/ (variance Ci) (variance S))) clustering))

   Returns a floating point density score of compactness - the smaller
   the better.
  "
  [distfn avgfn clustering]
  (let [S (apply set/union (xfold #(set (second %)) clustering))
        Svar (variance S :distfn distfn :avgfn avgfn :rfn xfold)
        Cvars (xfold (fn[[mi xis]]
                       (variance xis :distfn distfn :m mi))
                     1 clustering)]
    (mean (map #(/ % Svar) Cvars))))

(defn S-Dbw-index
  "Compute the 'SDbw' cluster validity index.  By several
   accounts (IEEE 2010 ICDM paper 'Understanding Internal Clustering
   Validity Measures', in particular) this is the most robust general
   internal validity measure across both data sets and clustering
   algorithms.  It is the default used by FIND-CLUSTERS.

   It accounts for both compactness of clusters and cluster
   separation.  It does this with a dual density measure:

   1. Scattering (see SCATT), which computes the average variance in
      clusters to the overall variance of the data set.

   2. Intercluster density (so called 'Dens_bw', see
      INTERCLUSTER-DENSITY), which computes an averaged density in the
      space between all cluster pairs.

   Note in particular that for convex sets, it is proved that the
   clustering which minimizes Scatt + Dens_bw, is the optimal
   clustering for the data set and algorithm pair (see Halkidi &
   Vazirgiannis, Clustering Validity Assesment: Finding Optimal
   partitioning of a Data Set, Proc ICDM 2001, pp187-194).

   DISTFN is the distance function for the data and AVGFN the 'mean'
   for the data.

   Returns a floating point number as score.  The smaller the better.
  "
  [distfn clustering & {:keys [avgfn] :or {avgfn mean}}]
   (+ (scatt distfn avgfn clustering)
      (intercluster-density distfn avgfn clustering)))



;;; (ns-unmap *ns* 'find-clusters)
(defmulti
  ^{:arglists
    '([coll & {:keys [distfn avgfn algo vindex]
             :or {algo kmeans++ vindex S-Dbw-index
                  distfn edist avgfn mean}}]
      [_ coll & {:keys [distfn avgfn algo vindex]
                 :or {algo kmeans++ vindex S-Dbw-index
                      distfn edist avgfn mean}}])}
  find-clusters
  "Computes the 'best' clustering of the data in collection COLL whose
   distances are given by DISTFN and means by AVGFN, as produced by
   ALGO and measured by VINDEX.  ALGO is the clustering algorithm,
   defaults to kmeans++, and VINDEX is the validity index measure,
   defaults to S-Dbw-index.

   'Best', here is as determined by vindex.  If the data has
   convex (natural/true) clusters, and is not skewed, the default
   kmeans++ with S-Dbw-index will return the optimal
   clustering (indeed, it will be or be extremely close to the
   natural/true clustering).

   If the data are not convex, kmeans++ is invalid (can only find
   convex clusters...).  If the data is heavily skewed, kmeans will
   not find the optimal (true) clusters (as it always tends to find
   'equal area' clusters.
  "
  (fn [& args] (first args)))


(defmethod find-clusters :default
  [coll & {:keys [distfn avgfn algo vindex]
           :or {algo kmeans++ vindex S-Dbw-index
                distfn edist avgfn mean}}]
  (let [data (set coll)]
    (loop [k 2
           prev-info [(double Integer/MAX_VALUE) []]]
      (let [[prev-score prev-clustering] prev-info
            clustering (algo k data :distfn distfn :avgfn avgfn)
            score (vindex distfn clustering :avgfn avgfn)]
        (if (> score prev-score)
          [(dec k) prev-score prev-clustering]
          (recur (inc k)
                 [score clustering]))))))


(defmethod find-clusters :global
  [_ coll & {:keys [distfn avgfn algo vindex]
             :or {algo kmeans++ vindex S-Dbw-index
                  distfn edist avgfn mean}}]
  (let [data (set coll)]
    (first
     (sort-by
      first <
      (for [k (range 2 (int (/ (count coll) 3)))
            :let [clustering (algo k data :distfn distfn :avgfn avgfn)
                  score (vindex distfn clustering :avgfn avgfn)]]
        [score clustering])))))


(comment
  (let [ents-seqs (->> "/home/jsa/Bio/FreqDicts/NewRFAM/RF00504-seed-NC.sto"
                       (#(read-seqs %1 :info :names))
                       (#(get-adjusted-seqs %1 700))
                       (take 30)
                       ((fn[[nm sq]] [(probs 7 sq) nm sq])) vec)
        ??? (combins 2 (range 0 30))]
    (clu/find-clusters :global xxx
                       :distfn ???
                       :avgfn ???))
  )


;;; ----------------- Ad Hoc Testing Stuff -------------------------------

(comment

(def data1 [2 3 5 6 10 11 100 101 102])
(def data (flatten (for [i (range 5)]
                     (let [j (int (rand 1000))
                           k (int (rand 10))]
                       (for [h (range k)]
                         (+ j (rand)))))))

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

(variance [[1 2 3] [3 2 1] [100 200 300] [300 200 100] [50 50 50]]
          :distfn vecdist :avgfn vecmean)

(std-deviation [[1 2 3] [3 2 1] [100 200 300] [300 200 100] [50 50 50]]
               :distfn vecdist :avgfn vecmean)

)