;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                                   W W W                                  ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011 Trustees of Boston College                            ;;
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

(ns edu.bc.bio.gaisr.www

  "Web server front end for Gaisr.  Based on an asyncrhonous evented
   http server (Aleph/netty), the Compojure route, file, and content
   type service frontend, and the Ring (monad like) server framework."

  (:require [clojure.contrib.string :as str]

            [ring.util.response :as ring-resp]
            [compojure.route :as route]
            [compojure.handler :as handler]
            [compojure.response :as response]

            [edu.bc.fs :as fs])

  (:use edu.bc.utils
        [edu.bc.log4clj :only [create-loggers log>]]
        edu.bc.bio.gaisr.actions

        clojure.contrib.math
        [clojure.pprint
         :only [cl-format]]

        [ring.middleware.cookies
         :only [wrap-cookies]]
        [ring.middleware.file-info
         :only [wrap-file-info]]
        [ring.middleware.stacktrace
         :only [wrap-stacktrace]]
        [ring.handler.dump
         :only [handle-dump]]
        [ring.middleware.etag
         :only [wrap-etag]]
        compojure.core
        lamina.core
        aleph.http)

  (:import (java.nio ByteBuffer CharBuffer)
           (java.io PushbackReader
                    InputStream
                    InputStreamReader
                    FileInputStream)))



;;; -------------------------------------------------------------------------;;;
;;; Web Server Stuff.  Uses Compojure with Aleph http evented server.


(defroutes main-routes
  (GET  "/" [] "<h1>Hello.  Server is up...</h1>")
  (GET  "/query" [& args] (ctx-query args))
  (GET  "/mlab/query" [:as reqmap & args] (ctx-query args reqmap))
  (GET  "/mlab/action" [& args] (action-request args))
  (POST "/mlab/action" [& args] (action-request args))
  (POST "/mlab/upld"   [:as reqmap & args] (upload-file args reqmap))
  (GET  "/mlab/features" [& args] (get-features args))

  (GET ["/:filename" :filename #".*"] [filename]
       (ring-resp/file-response filename {:root "./resources"}))
  (route/resources "/")

  ;;(route/not-found "No Such Page")
  (ANY "*" {:keys [uri params form-params] :as request}
       (fall-through-catcher uri params form-params request)))




(defn wrap-status [handler]
  (fn [request]
    (let [response (handler request)]
      (assoc response :status
             (if-let [err (:error response)]
               err
               200)))))


(defn get-content-type [response]
  (-> response :headers
      ((fn [headers]
         (or (get headers "Content-Type")
             (get headers "content-type"))))))

(defn wrap-content-type [handler]
  (fn [request]
    ;;(log> "rootLogger" :info "WCT: ~S" request)
    (let [request
          (if (= (:request-method request) :get)
            request
            (let [content-type (get-in request [:headers "content-type"])
                  charset (re-find #"charset=.*$" content-type)
                  character-encoding (if charset
                                       (second (str/split #"=" charset))
                                       "UTF-8")
                  content-length (get-in request [:headers "content-length"])
                  request (assoc request
                            :content-length (Integer. content-length)
                            :content-type content-type
                            :character-encoding character-encoding)]
              request))
          response (handler request)]
      (if (get-content-type response)
        response
        (-> response
            (assoc-in [:headers "content-type"]
                      (if (= (:type response) :json)
                        "application/json"
                        "text/html"))
            (dissoc :type))))))




(defn wrap-cache-control [handler]
  (fn [request]
    (let [response (handler request)]
      (assoc-in response
                [:headers "Cache-Control"]
                ["no-transform" "no-cache"]))))


(def req-log (atom ""))
(def req-body-log (atom ""))
(def log-req? (atom false))
(def log-resp? (atom false))

(defn log-request [handler]
  (fn [request]
    (when @log-req?
      (swap! req-log (fn[_ rq] rq) request)
      (swap! req-body-log (fn[_ b] b) (:body request))
      (log> "rootLogger" :info "REQ: ~S" request))
    (handler request)))

(defn log-response [handler]
  (fn [request]
    (let [response (handler request)]
      (when @log-resp?
        (log> "rootLogger" :info "~A" response))
      response)))


;;; The main app is a pipeline of middleware wrappers to main site
;;; ports main-routes which dispatches to main processors (ctx-query,
;;; feature-query, blast, etc)
;;;
(def www-gaisr-ring
     (-> (handler/site main-routes)
         log-request
         ;;wrap-stacktrace
         wrap-cookies
         wrap-status
         wrap-file-info
         wrap-content-type
         wrap-cache-control
         wrap-etag
         log-response
         ))


(defn www-gaisr [ch request]
  ;;(if (nil? (:body request)) ; remove when new Aleph available...
  ;;  (enqueue ch (www-gaisr-ring request))
  (run-pipeline
   (request-body->input-stream request)
   www-gaisr-ring
   #(enqueue ch %)))


;;; Basic Aleph server app scaffolding.  The main difference between
;;; this and normal Ring is simply that we need a _channel_ as well as
;;; the request since this is going to work _asynchronously_.  The
;;; part that is the "normal" Ring app is body request handler.
;;;
(defn gaisr-main [channel request]
  (www-gaisr channel request))


(def stop-server-closure
     (atom (fn[] (raise :type :no-server :args "Server not yet started!"))))

(defn stop-gaisr-server []
  (@stop-server-closure))

(defn start-gaisr-server [& {port :port :or {port 8080}}]
  (swap! stop-server-closure
         (fn[_] (start-http-server gaisr-main {:port port}))))

;;;(stop-gaisr-server)
;;;(start-gaisr-server)
