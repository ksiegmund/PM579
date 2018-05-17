## intro.R

getwd()
# for people running this on lab machines
#setwd("U:")
#getwd()

2+5
x=2+5
x
x=c(1,3,5,7)
x
y=c(2,4,6,8)
x+y
z=x+y
z

# function calls are followed by ()
#quit()

y=c(x,0,x)
y
2*x
2*x+y+1

x <- array(1:20, dim=c(4,5)) # Generate a 4 by 5 array (or matrix)

Lst <- list(name="Fred", wife="Mary", no.children=3,
child.ages=c(4,7,9))
Lst
Lst$child.ages

make=c("Honda","Toyota","Ford","GM")
year=c(1997,2002,2007,2005)
cars=data.frame(make,year)
cars
cars[1,]
cars$make
cars$year

##
##  Gene Expression
##
# Import C42b cell line data from GEO
#Illumina HumanHT-12 V4.0 expression beadchip
#GPL10558; 47231 rows
source("http://bioconductor.org/biocLite.R")
if(!require("GEOquery")) biocLite("GEOquery")
library(GEOquery)
# The following line will download from GEO. However, this will
# be extremely slow if we all do it simultaneously over the wifi.
# Therefore I ran the command and saved the object as an R data set
# that we can load with the following command.
#gse31873=getGEO('GSE31873',GSEMatrix=TRUE)
load("data/JBC 2012/gse31873.rda")
show(gse31873)
c42b=gse31873$GSE31873_series_matrix.txt.gz
ec42b=exprs(c42b)    # exprs() accesses the gene expression values
ec42b[1:4,1:3] 
names(pData(c42b))   # pData() accesses the sample annotation information
# Check if the sample annotation data and gene expression data
# are ordered identically
identical(as.character(pData(c42b)$geo_accession),
          colnames(ec42b))
pData(c42b)$title    # sample treatment information

##
## MDS plot
##
# use treatment to label points on the plot
trtmat<- matrix(
     unlist(strsplit(as.character(pData(c42b)$title),split=" ")),
            nrow=24,byrow=TRUE)
trtmat<- data.frame(trtmat)
colnames(trtmat)<-c("cellline","TrtTime","rx","Rep")
library(limma)
plotMDS(ec42b,labels=paste(trtmat$TrtTime, 
                            unclass(trtmat$Rep), sep="_"),
        col=unclass(trtmat$TrtTime),
        xlim = c(-1.5,1.5), ylim=c(-1,1),
        main="MDS plot") #color by type

##
## Heatmap
##
#source("http://bioconductor.org/biocLite.R")
if(!require("matlab")) biocLite("matlab")
if(!require("tidyverse")) install.packages("tidyverse")
library(matlab)   # this library let's us use blue-red color spectrum for heatmap  (jet.colors)
library(tidyverse)      # this library has the recode command

## assign colors to sample annotation variables
clab<- matrix(unlist(strsplit(as.character(trtmat$TrtTime),split="_")),
               ncol=2,byrow=T)
colnames(clab)<-c("treatments","hour")
level_key <- list(siNS = "pink", sip300 = "orange", siCBP = "purple", '0hr' = "grey", '16hr' = "lightgreen")
clab[,1]<- recode(clab[,1], !!!level_key)
clab[,2]<- recode(clab[,2], !!!level_key)
clab

top500var=which(rank(apply(ec42b,1,IQR))>nrow(ec42b)-500)
length(top500var)

source("Intro to R/heatmap.3.R")
par(cex.main=1)
heatmap.3(ec42b[top500var,], scale="none",  margins=c(2,10), Rowv=TRUE, Colv=TRUE,
          ColSideColors=clab,NumColSideColors=4, labCol=FALSE, labRow=NA,cexCol=1.8,
          symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", KeyValueName="log2(Expr)",
          col=jet.colors(32))
legend("topright",
       c("NS","CBP","p300","0hr","16hr"),cex=1.2,
       fill=c("pink","purple","orange","grey","lightgreen"), border=FALSE,bty="n")


sessionInfo()
