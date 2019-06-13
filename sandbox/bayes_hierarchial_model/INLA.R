library(INLA)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)


Hosts <- read.csv("~/Downloads/CC-INLA-master/HostCaptures.csv", header = T)
head(Hosts)

substr(names(Hosts), 1, 1) <- toupper(substr(names(Hosts), 1, 1)) # Giving the host names capital letters

phen <- c("Grid", "ID", "Easting", "Northing") # Base columns with spatial information we'll need

resp <- "Parasite.count" # Response variable

covar <- c("Month", # Julian month of sampling
           "Sex", # Sex
           "Smi", # Body condition
           "Supp.corrected", # Nutrition supplementation
           "Treated") # Treatment

TestHosts <- na.omit(Hosts[, c(phen, resp, covar)]) # Getting rid of NA's, picking adults
# We are using the [] to subset and only extract specific columns

# Turning variables into factors
TestHosts$Month <- as.factor(TestHosts$Month)
TestHosts$Grid <- as.factor(TestHosts$Grid)

TestHosts$Parasite.count <- round(TestHosts$Parasite.count) # Parasite counts should be integers

table(table(TestHosts$ID)) # Enough repeat samples for a mixed model?



# Setting up a custom theme
THEME <- theme(axis.text.x = element_text(size = 12,colour = "black"),
               axis.text.y = element_text(size = 12, colour = "black"),
               axis.title.x = element_text(vjust = -0.35),
               axis.title.y = element_text(vjust = 1.2)) + theme_bw()

samp_locations <- ggplot(TestHosts, aes(Easting, Northing)) + 
    geom_jitter(aes(colour = factor(Grid))) + coord_fixed() + 
    THEME + 
    labs(colour = "Grid")


# Specify the formula
f0.1 <- as.formula(paste0(resp, " ~ ", # Response first
                          paste(covar, collapse = " + ") # Collapse the vector of covariates
))

# Run the model
IM0.1  <- inla(Parasite.count ~ Month + Sex + Smi + Supp.corrected + Treated, 
               family = "nbinomial", # Specify the family. Can be a wide range (see r-inla.org).
               data = TestHosts) # Specify the data

# Run the model # (This is the same thing)
IM0.1  <- inla(f0.1, 
               family = "nbinomial", # Specify the family. Can be a wide range (see r-inla.org).
               data = TestHosts) # Specify the data

# Then with an ID random effect ####

f0.2 <- as.formula(paste0(resp, " ~ ", 
                          paste(covar, collapse = " + "), 
                          " +  f(ID, model = 'iid')")) # This is how you include  a typical random effect.

IM0.2  <- inla(f0.2, 
               family = "nbinomial",
               data = TestHosts) 

summary(IM0.1)
summary(IM0.2)

dat <- data.frame(y=runif(100,0,1), x=rbinom(100, 1, 0.5))

INLA::inla(y ~ x, data = dat, verbose = T)

