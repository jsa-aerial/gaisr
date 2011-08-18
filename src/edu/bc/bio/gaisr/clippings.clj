
;;; ***NOTE: No longer needed or used.  Now fixed and handled by
;;; ***Aleph!!
;;;
;;; Low level "nasty stuff" to handle form bodies coming from netty,
;;; since Ring does not have an out of the box handler for netty yet.
;;; Maybe this should be fixed in Aleph???
;;;
;;; Netty form bodies come over as _arrays_(!?!?) of ByteBuffers, but
;;; (to make it even odder) with _only one element_!  Ugh.  Ring
;;; expects an input stream for the form body.  So, we place a handler
;;; ahead of all others which transforms the Netty body into an input
;;; stream...
;;;
(defn buffer-stream
  "Returns an InputStream for a ByteBuffer."
  [#^ByteBuffer buf]
  (proxy [InputStream] []
    (available [] (.remaining buf))
    (read
      ([] (if (.hasRemaining buf) (.get buf) -1))
      ([dst offset len] (let [actlen (min (.remaining buf) len)]
                          (.get buf dst offset actlen)
                          (if (< actlen 1) -1 actlen))))))

(defn handle-aleph-body [request]
  (if (not (:body request))
    request
    (let [content-type (get-in request [:headers "content-type"])
          charset (re-find #"charset=.*$" content-type)
          character-encoding (if charset
                               (second (str/split #"=" charset))
                               "UTF-8")
          body (:body request)
          body (if (channel? body)
                 (do
                   (receive-all body
                     (fn[chunk]
                       (println chunk)
                       (set! req-body (conj req-body chunk))))
                   (buffer-stream (first (reverse req-body))))
                 (buffer-stream (first body)))
          request (assoc request
                    :content-type content-type
                    :character-encoding character-encoding
                    :body body)]
      request)))


