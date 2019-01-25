library(dagitty)
library(ggdag)

confounder_triangle(x = "Exposure", 
                    y = "Outcome", 
                    z = "Counfounder") %>% 
  ggdag_dconnected(text = FALSE, use_labels = "label")

pred <- dagitty( "dag {
               X1 -> Y
               X2 -> Y
               X3 -> Y
               X4 -> Y
               X5 -> Y
               X6 -> Y
               }")


pred <- dagify(    y ~ x1,
                   y ~ x2,
                   y ~ x3,
                   y ~ x4,
                   y ~ x5,
                   y ~ x6,
                  exposure = "x1",
                  outcome = "y")
tidy_dagitty(pred)
ggdag::ggdag_adjustment_set(pred) + 
  theme_dag_blank() + 
  geom_dag_text(color="black") + 
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
        )


