(in-ns 'edu.bc.bio.canned-jobs)

(def *jid* (add *jobs*
                '[(t1-blast "/data2/BioData/BlastDBs/fasta-test.fna"
                            "/data2/BioData/fasta-test.fna.blast"
                            "/data2/BioData/BlastDBs/other_genomic.00")
                  (t2-blast :t1)]
                jfn-blast))