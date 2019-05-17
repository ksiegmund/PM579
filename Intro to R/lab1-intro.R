## intro.R

getwd()

2+5
x <- 2+5
x
x <- c(1,3,5,7)
x
y <- c(2,4,6,8)
x+y
z <- x+y
z

# function calls are followed by ()
#quit()

y <- c(x,0,x)
y
2*x
2*x+y+1

x <- array(1:20, dim=c(4,5)) # Generate a 4 by 5 array (or matrix)

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

##
##  Gene Expression
##
# Import C42b cell line data from GEO
#Illumina HumanHT-12 V4.0 expression beadchip
#GPL10558; 47231 rows
if(!require("GEOquery")){BiocManager::install("GEOquery")}
library(GEOquery)
# The following line will download from GEO. However, this will
# be extremely slow if we all do it simultaneously over the wifi.
# Therefore I ran the command and saved the object as an R data set
# that we can load with the following command.
#gse31873=getGEO('GSE31873',GSEMatrix=TRUE)
load("data/JBC 2012/gse31873.rda")

show(gse31873)
gse <- gse31873$GSE31873_series_matrix.txt.gz
# I'm going to create a list object, jbcdat, that consolidates
# 1. expression measures   exprs()
# 2. sample annotation     pData()    p is short for phenotype
# 3. feature annotation    fData()
jbcdat <- NULL
jbcdat$E <- exprs(gse)    # exprs() accesses the gene expression values
jbcdat$E[1:4,1:3] 
names(pData(gse))   # pData() accesses the sample annotation information
names(fData(gse))   # fData() accesses the feature annotation information
# Check if the sample annotation data and gene expression data
# are ordered identically
identical(as.character(pData(gse)$geo_accession),
                colnames(jbcdat$E))
pData(gse)$title        # sample treatment information
pData(gse)$description  # array name

##  create a sample annotation dataframe named targets
trtvars <- matrix(
            unlist(strsplit(as.character(pData(gse)$title),split=" ")),
            nrow=24,byrow=TRUE)
trtvars <- data.frame(trtvars)
colnames(trtvars) <- c("cellline","trt_time","rx","rep")
head(trtvars)
targets <- matrix(
              unlist(strsplit(as.character(trtvars$trt_time),split="_")),
              nrow=24,byrow=TRUE)
targets <- data.frame(targets)
colnames(targets) <- c("treatments","hour")
rownames(targets) <- pData(gse)$description
head(targets)
jbcdat$targets <- data.frame(targets,rep = trtvars$rep, type = trtvars$trt_time, 
                       trt_time_rep = paste(trtvars$trt_time,trtvars$rep,sep="_"))
head(jbcdat$targets)
names(jbcdat)
rm(targets)

# Now let's save some array feature annotation information 
jbcdat$genes <- data.frame(fData(gse)[,c("ID","Entrez_Gene_ID","Symbol","Chromosome")])
head(jbcdat$genes,n=12)
names(jbcdat)

# Let's save this data set for future use
jbcdir=c("data/JBC 2012")
save(jbcdat,file = file.path(jbcdir,"jbcdat.rda"))


##
## MDS plot
##
# use treatment to label points on the plot
library(limma)

limma::plotMDS(jbcdat$E,labels=jbcdat$targets$trt_time_rep,
        col=unclass(jbcdat$targets$type),
        xlim = c(-1.5,1.5), ylim=c(-1,1),
        main="MDS plot") #color by type

##
## Heatmap
##
if(!require("ComplexHeatmap")){BiocManager::install("ComplexHeatmap")}
if(!require("matlab")){BiocManager::install("matlab")}
library(ComplexHeatmap)
library(matlab)   # this library let's us use blue-red color spectrum for heatmap  (jet.colors)

# get the row (feature) number for the 500 features with largest IQR
fiqr <- apply(jbcdat$E,1,IQR)
fsum <- data.frame(fiqr, rfiqr = rank(-fiqr))
head(fsum)
fidx <- which( fsum$rfiqr <= 500)
length(fidx)

# column heatmap annotation
colha = HeatmapAnnotation(df = jbcdat$targets[,c("treatments","hour")],
                          col = list(treatments = c(siNS = "pink", 
                                                   siCBP = "purple",
                                                   sip300 = "orange"),
                                           hour = c('0hr' = "grey",
                                                    '16hr' = "lightgreen")
                                     ), 
                          which = "column")
# heatmap
ht = Heatmap(jbcdat$E[fidx,], column_title = "Samples",
              row_title = "Features", 
              name = "log2(Expr)",
              col = jet.colors(32), 
              top_annotation = colha,
              show_column_names = FALSE,
              show_row_names = FALSE)
draw(ht)

sessionInfo()
