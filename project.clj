(defproject net.clojars.oyvinht/markov-chains "0.9.3"
  :description "A library for working with Markov Chains."
  :url "https://github.com/oyvinht/markov-chains"
  :license {:name "EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0"
            :url "https://www.eclipse.org/legal/epl-2.0/"}
  :dependencies [[com.clojure-goes-fast/clj-async-profiler "0.5.0"]
                 [criterium "0.4.6"]
                 [org.clojure/clojure "1.10.1"]]
  :jvm-opts ["-Djdk.attach.allowAttachSelf" "-XX:+UnlockDiagnosticVMOptions" "-XX:+DebugNonSafepoints"]
  :repl-options {:init-ns markov-chains.core}
  :profiles {:uberjar {:aot :all}}
  :plugins [[lein-nodisassemble "0.1.3"]]
  )
