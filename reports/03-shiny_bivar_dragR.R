library(tidyverse)
library(dragulaR)
library(shiny)

#........................
# setup functions
#........................
rotate <- function(x){
  if(length(x) == 4){
    ret <- t(apply(x, 2, rev))
    return(ret)
  } else if(length(x) == 8){
    for(i in 1:(log2(length(x))-1)){
      x[,,i] <- t(apply(x[,,i], 2, rev)) # note, oddly the row/col names don't get dragged along in the array. oh well
    }
    return(x)
  } else{
    return(x) # do nothing because we are now in too many dimensions for this to work simply
  }
}

shiny_print_epi2by2 <- function(x){
  if(x$method == "case.control" & x$n.strata == 1){ ## case.control --- single strata
    
    p1 <- with(x$res, {
      
      sprintf("\nOdds ratio (Wald)            %.2f (%.2f, %.2f)",
              OR.strata.wald[[1]], 
              OR.strata.wald[[2]],
              OR.strata.wald[[3]]
      )})
    p2 <- with(x$res, {
      sprintf("\nCrude Chi-sq (stat, df, p)   %.2f, %.2f, %.2e",
              chisq.crude[1],
              chisq.crude[2],
              chisq.crude[3]
      )})
    ret <- c(p1,p2)
    return(ret)
    
  } # close out if == nstrata=1 
  else if(x$method == "case.control" & x$n.strata > 1){
    
    p1 <- with(x$res, {
      
      sprintf("\nOdds ratio (crude)             %.2f (%.2f, %.2f)",
              OR.crude.wald[1],
              OR.crude.wald[2],
              OR.crude.wald[3]
      )})
    p2 <- with(x$res, {
      sprintf("\nOdds ratio (M-H)               %.2f (%.2f, %.2f)",
              OR.mh.wald[1],
              OR.mh.wald[2],
              OR.mh.wald[3]
      )})
    p3 <- with(x$res, {
      sprintf("\nRatio of OR (crude:M-H)        %.2f",
              round(OR.crude.wald[1] / OR.mh.wald[1], digits = 2)
      )})
    p4 <- with(x$res, {
      sprintf("\nCrude Chi-sq (stat, df, p)     %.2f, %.2f, %.2e",
              chisq.crude[1],
              chisq.crude[2],
              chisq.crude[3]
      )})
    p5 <- with(x$res, {
      sprintf("\nStrat Chi-sq (stat, df, p)     %.2f, %.2f, %.2e",
              chisq.strata[[1]],
              chisq.strata[[2]],
              chisq.strata[[3]]
      )})
    ret <- c(p1,p2,p3,p4,p5)
    return(ret)
    
  } # close out ifelse for n.strata>1
  else {
    stop("This simple function does not support that input type")
  }
} # end function      

#........................
# import data and covars
#........................
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
dt$adm1name_fctm <- factor(dt$adm1name)
dt$hv001_fctm <- factor(dt$hv001)
covars <- colnames(dt)[grepl("_fctm|_fctb|_cont", colnames(dt))]

makeElement <- function(data, name)
{
  div(style = "border-width:1px;border-style:solid;",
      drag = name,
      div(class = "active title", name),
      div(class = "active content", p(sprintf("Class: %s", class(data[[name]])))))
}

ui <- fluidPage(
  titlePanel("Drag and drop covariates"),
  
  fluidRow(style = "margin: 15px;",
           column(3,
                  h3("Drag from here:"),
                  div(id = "Available", style = "min-height: 2400px;",
                      lapply(covars, makeElement, data = dt))
           ),
           column(3,
                  h3("Drop here:"),
                  div(id = "Model", style = "min-height: 2400px;")
           ),
           column(6,
                  plotOutput("plot"),
                  tableOutput("tableone_kable"),
                  br(),
                  tableOutput("table"),
                  br(),
                  tableOutput("prev"),
                  verbatimTextOutput("print")
           )
  ),
  dragulaOutput("dragula")
  
)

server <- function(input, output) {
  
  output$dragula <- renderDragula({
    dragula(c("Available", "Model"))
  })
  
  output$plot <- renderPlot({
    req(input$dragula)
    state <- dragulaValue(input$dragula)
    validate(need(length(state$Model) > 1, message = " "))
    
    plot(dt[,state$Model])
  })
  
  output$tableone_kable <- function(){
    req(input$dragula)
    state <- dragulaValue(input$dragula)
    validate(need(length(state$Model) > 1, message = "Please select at least two variables."))
    tableone::CreateTableOne(vars = state$Model[2:length(state$Model)],
                               strata = state$Model[1],
                               data = dt,
                               includeNA = T) %>%
      tableone::kableone(., "html") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = T) %>% 
      kableExtra::add_header_above(c("Table One Statistics" = 5))
    
  }
  
  
  output$table <- function(){
    state <- dragulaValue(input$dragula)
    validate(need(length(state$Model) > 1, message = " "))
    as.data.frame(xtabs(paste0("~", paste(state$Model, collapse = "+")), data=dt)) %>% 
      tableone::kableone(., "html") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = T) %>% 
      kableExtra::add_header_above(c("Covariate Counts" = 3))
  }
  
  output$prev <- function(){
    state <- dragulaValue(input$dragula)
    validate(need(length(state$Model) > 1, message = " "))
    ret <- epiR::epi.2by2(dat = rotate(rotate( 
      xtabs(paste0("~", paste(state$Model, collapse = "+")), data=dt)
      )),
                          method = "case.control",
                          conf.level = 0.95, homogeneity = "breslow.day")
    
    ret <- as.data.frame(ret$tab)
    ret %>% 
      tableone::kableone(., "html") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = T) %>% 
      kableExtra::add_header_above(c("Two-by-Two & Summary Effect Estimates" = 6))
    
  }
  
  output$print <- renderText({
    state <- dragulaValue(input$dragula)
    validate(need(length(state$Model) > 1, message = " "))
    ret <- epiR::epi.2by2(dat = rotate(rotate(
      xtabs(paste0("~", paste(state$Model, collapse = "+")), data=dt)
      )), # annoying epiR setup
                         method = "case.control",
                         conf.level = 0.95, homogeneity = "breslow.day")
    ret <- shiny_print_epi2by2(ret)
    paste(ret)
    
  })
  
}

shinyApp(ui = ui, server = server)

