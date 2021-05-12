library(dagitty)

# m dag
dag <- dagitty::downloadGraph(x = "dagitty.net/m8uo5fS")
adj <- dagitty::adjustmentSets(dag, 
                               exposure = "Exposure",
                               outcome = "Outcome",
                               type = "canonical",
                               effect = "total")
adj

# m-plus dag
dag <- dagitty::downloadGraph(x = "dagitty.net/mnAWT6k")
adj <- dagitty::adjustmentSets(dag, 
                               exposure = "Exposure",
                               outcome = "Outcome",
                               type = "canonical",
                               effect = "total")
adj