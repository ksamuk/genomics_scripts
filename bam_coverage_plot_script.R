source("http://bioconductor.org/biocLite.R")
biocLite("chipseq")

library("Rsamtools")
library("IRanges")
library("chipseq")
bam <- scanBam("whtstbk_gbs_2015_L1_brds_GC14.realign.bam")

bam <- scanBam("whtstbk_gbs_2015_L1_brds_GC24.realign.bam")

names(bam[[1]])
table(bam[[1]]$flag)

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) 
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}

.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

test <- subset(bam_df, rname == 'groupIV')
test <- subset(bam_df, rname == 'groupXV')

no.nas.tmp <- data.frame(pos = test$pos, qwidth = test$qwidth)
no.nas.tmp <- no.nas.tmp[complete.cases(no.nas.tmp), ]

dat.tmp<-no.nas.tmp

bam.ir  <- with(dat.tmp, IRanges(start = pos, width = qwidth))
plotRanges(bam.ir)
intervals<-ranges(coverage(bam.ir))
olaps<- findOverlaps(narrow(intervals,width(intervals)),bam.ir)

bam.ir[tabulate(queryHits(olaps), queryLength(olaps))>25000,]

bins <- disjointBins(IRanges(start(bam.ir ), end(bam.ir ) + 1))
# [1] 1 2 1 2 3 1 1

dat <- cbind(as.data.frame(bam.ir ), bin = bins)

library(ggplot2)
ggplot(dat) + 
  geom_rect(aes(xmin = start, xmax = end,
                ymin = bin, ymax = bin + 0.9)) +
  theme_bw()

autoplot(Cov)


#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)

#use chr22 as an example
#how many entries on the negative strand of chr22?
table(bam_df$rname == 'groupIV' & bam_df$flag == 16)
# FALSE    TRUE 
#3875997   24413

#function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

#test neg function with subset of chr22
test <- subset(bam_df, rname == 'groupIV')
dim(test)
#[1] 56426    13
table(apply(as.data.frame(test$flag), 1, check_neg))
#number same as above
#FALSE  TRUE 
#32013 24413

#function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[3] != 1 & intToBits(x)[1] == 0){
    return(T)
  } else {
    return(F)
  }
}

#check pos function
table(apply(as.data.frame(test$flag), 1, check_pos))
#looks OK
#FALSE  TRUE 
#24413 32013

#store the mapped positions on the plus and minus strands
chr22_neg <- bam_df[bam_df$rname == 'groupIV' &
                      apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos'
                    ]
length(chr22_neg)
#[1] 24413
chr22_pos <- bam_df[bam_df$rname == 'groupIV' &
                      apply(as.data.frame(bam_df$flag), 1, check_pos),
                    'pos'
                    ]
length(chr22_pos)
#[1] 32013

#calculate the densities
chr22_neg_density <- density(chr22_neg)
chr22_pos_density <- density(chr22_pos)

#display the negative strand with negative values
chr22_neg_density$y <- chr22_neg_density$y * -1

plot(chr22_neg_density,
     main = "Coverage plot of mapped CAGE reads",
     xlab = "Chromosome IV",
     col = 'blue',
     lwd=2.5)
lines(chr22_neg_density, lwd=2.5, col = 'red')