# pca of filtered vcf

source("http://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")

install.packages("pegas")
library("pegas")
library("adegenet")
library("VariantAnnotation")


# read in vcf file
plink.file <- list.files(pattern=".raw")
plink.map <- list.files(pattern=".map")

gbs2015 <- read.PLINK(plink.file, map.file=plink.map)

#get pop ids

ids.raw <- as.character(pop(gbs2015))

strip_prefix <- function(x){
  split <- strsplit(x, split="_")
  return(split[[1]][length(split[[1]])])
}

ids <- sapply(ids.raw, strip_prefix)
ids <- unname(ids)
pop.ids <- gsub("\\d", "", ids)
pop.ids <- gsub("\\W", "", pop.ids)
pop(gbs2015) <- as.factor(pop.ids)

is.genlight(gbs2015)

gbs.df <- as.data.frame(gbs2015)

gbs.genind <- df2genid(gbs.df)

# pca plot

# 3 and 2 retained? 

gbs.dapc <- dapc(gbs2015)

#3 axes
gbs.pca <- glPca(gbs2015)

scatter(gbs.pca, xax=1, yax=2,grid=FALSE, label=ids, ratio=0.6, posi="none", clabel=0.5, col="red")

