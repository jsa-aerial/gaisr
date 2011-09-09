(defproject gaisr "1.0.0-SNAPSHOT"
  :description "Genomic Analysis Illuminating Stuctured RNA pipeline, scoring, curation system for ncRNA"
  :dependencies [[org.clojure/clojure "1.2.1"]
                 [org.clojure/clojure-contrib "1.2.0"]
                 [clojure-csv "1.3.0"]
                 [clj-shell "0.1.0"]

                 ;;[clj-native "0.9.1-SNAPSHOT"]
                 [org.clojars.nathell/clojure-jna "1.0.0-SNAPSHOT"]

                 [ring "0.3.5"]
                 [compojure "0.6.0-RC4"]
                 [aleph, "0.2.0-beta1"]

                 [org.bituf/clj-dbcp "0.2"]
                 [c3p0/c3p0 "0.9.1.2"]
                 [mysql/mysql-connector-java "5.1.15"]
                 [commons-net "2.2"]

                 [org.clojure.contrib/macro-utils "1.3.0-alpha4"]
                 [core.logic "0.6.1-SNAPSHOT"]

                 [swingrepl "1.0.0-SNAPSHOT"]
                 [incanter/incanter-core "1.2.3"]
                 [incanter/incanter-io "1.2.3"]
                 [incanter/incanter-charts "1.2.3"]
                 [incanter/incanter-processing "1.2.3"]
                 [incanter/incanter-mongodb "1.2.3"]
                 [incanter/incanter-pdf "1.2.3"]
                 [incanter/incanter-latex "1.2.3"]
                 [incanter/incanter-excel "1.2.3"]
                 ]

  :dev-dependencies [[swank-clojure "1.3.2"]
                     [lein-clojars "0.6.0"]
                     [jline "0.9.94"]
                     ;;[mycroft/mycroft "0.0.2"]
                     ]

  :repositories {"commons-releases" "http://bizdirusa.com/mirrors/apache"
                 "c3p0-releases"    "http://mirrors.ibiblio.org/pub/mirrors/maven2"})
