## intro.R

getwd()

# taken from  http://cran.r-project.org/doc/manuals/R-intro.pdf

# Expression 

2+5

# Assignment

x <- 2+5
x
x <- c(1,3,5,7)
x
y <- c(2,4,6,8)

# Vector arithmetic

x+y
z <- x+y
z

y <- c(x,0,x)
y
2*x
2*x+y+1

# Arrays

x <- array(1:20, dim=c(4,5)) # Generate a 4 by 5 array (or matrix)

# Lists

Lst <- list(name="Fred", wife="Mary", no.children=3,
            child.ages=c(4,7,9))
Lst
Lst$child.ages

make <- c("Honda","Toyota","Ford","GM")
year <- c(1997,2002,2007,2005)
cars <- data.frame(make,year)
cars
cars[1,]
cars$make
cars$year

cars
is.character(cars$make)
is.factor(cars$make)
as.factor(cars$make)
levels(as.factor(cars$make))
as.numeric(as.factor(cars$make))


# Modern data frame

if(!require("tidyverse")) {install.packages("tidyverse")}
library("tidyverse")
cars <- tibble(make,year)

# Load a saved Gene expression object
jbcdir=c("data/JBC 2012")
load(file = file.path(jbcdir,"jbcdat.rda"))

names(jbcdat)
dim(jbcdat$E)
dim(jbcdat$targets)
dim(jbcdat$genes)

jbcdat$E[1:4,1:4]
head(jbcdat$targets)
head(jbcdat$genes)

identical(rownames(jbcdat$targets),colnames(jbcdat$E))
identical(rownames(jbcdat$E),rownames(jbcdat$genes))

# Dimension reduction technique
# Multi-dimensional scaling plot
if(!require("limma")) {BiocManager::install("limma")}

library(limma)
limma::plotMDS(jbcdat$E,pch=16,
               col=unclass(jbcdat$targets$type),
               xlim = c(-1.5,1.5), ylim=c(-1,1),
               main="MDS plot") #color by type
legend(-1.5,1,levels(jbcdat$targets$type),
        pch=16,col=order(levels(jbcdat$targets$type)),
       cex=.75)

sessionInfo()
