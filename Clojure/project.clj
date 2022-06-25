(defproject net.clojars.oyvinht/markov-chains "0.9.3"
  :description "A library for working with Markov Chains."
  :url "https://github.com/oyvinht/markov-chains"
  :license {:name "MIT"
            :url "LICENCE.TXT"}
  :dependencies [[com.clojure-goes-fast/clj-async-profiler "0.5.0"]
                 [criterium "0.4.6"]
                 [org.clojure/clojure "1.10.1"]]
  :jvm-opts ["-Djdk.attach.allowAttachSelf" "-XX:+UnlockDiagnosticVMOptions" "-XX:+DebugNonSafepoints"]
  :repl-options {:init-ns markov-chains.core}
  :profiles {:uberjar {:aot :all}}
  :plugins [[lein-nodisassemble "0.1.3"]]
  )
