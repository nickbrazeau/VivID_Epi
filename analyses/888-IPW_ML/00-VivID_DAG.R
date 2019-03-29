library(tidyverse)
library(dagitty)
library(ggdag)


dag <- dagitty::downloadGraph(x = "dagitty.net/mjKNQhB")
  
# DEFENSE OF DAG


ggdag::ggdag_adjustment_set(dag) + 
  theme_dag_blank() + 
  geom_dag_text(color="black") + 
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
        )


