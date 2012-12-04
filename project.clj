(defproject gaisr "1.2.0-SNAPSHOT"
  :description "GAISR - Genomic Analysis Illuminating Structured RNA"

  :dependencies [[org.clojure/clojure "1.5.0-beta1"]
                 [mlab.jars/clojure-contrib "1.2.0-mlab"]

                 [log4j/log4j "1.2.13"] ; base java log4j

                 [robert/bruce "0.7.1"] ; retries with decays, etc.
                 [slingshot "0.10.3"]   ; enhanced exceptions
                 [mlab.jars/clojure-csv "1.3.0-mlab"]
                 [clj-shell "0.1.0"]

                 [org.clojure.contrib/macro-utils "1.3.0-alpha4"]

                 ;;[clj-native "0.9.1-SNAPSHOT"]
                 [net.java.dev.jna/jna "3.4.0"]
                 [org.clojars.nathell/clojure-jna "1.0.0-SNAPSHOT"]

                 [ring "1.1.6"]
                 [org.clojars.jws/ring-etag-middleware "0.1.2-SNAPSHOT"]
                 [compojure "1.1.3"]
                 [aleph, "0.2.2"]

                 [org.clojure/java.jdbc "0.2.3"]
                 ;;[org.bituf/clj-dbcp "0.2"]
                 [commons-dbcp/commons-dbcp "1.4"] ; need this with mlab dbcp...
                 [mlab.jars/clj-dbcp "0.2-mlab"]
                 ;;[clj-dbcp "0.8.0"] ; Way too many changes (and dumber!)
                 [c3p0/c3p0 "0.9.1.2"]
                 [mysql/mysql-connector-java "5.1.15"]

                 ;;[clj-http "0.2.0"]
                 [org.clojars.ghoseb/enlive "1.2.0-alpha1"]
                 [commons-net "3.1"]

                 [swingrepl "1.3.0"]
                 [incanter "1.4.0"]

                 ;; SVM stuff
                 [nz.ac.waikato.cms.weka/weka-stable "3.6.8"]
                 [tw.edu.ntu.csie/libsvm "3.1"]
                 ]

  :dev-dependencies [[jline "0.9.94"]
                     [mlab.jars/swank-clojure "1.5.0-sd-mlab-col-hack"]
                     ]

  :repositories
  {"commons-releases" "http://bizdirusa.com/mirrors/apache"
   "c3p0-releases"    "http://mirrors.ibiblio.org/pub/mirrors/maven2"})
