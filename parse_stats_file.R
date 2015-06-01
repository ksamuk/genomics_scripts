# deep parsing of stats file

library(dplyr)
library(ggplot2)
library(data.table)
library("magrittr")

file.name <- "whtcmn_group1_slidingwindow_75000.txt"
file.name.snps <- "whtcmn_group1_stats.txt"
wht.stats <- read.table(file.name, header=TRUE)
wht.stats.snp <- fread(file.name.snps, header=TRUE)

# windows 
# set negative FSTs to NA
# this is sketchy, but for now i'm doing it
wht.stats[!is.na(wht.stats$Fst) & wht.stats$Fst<0, ]$Fst <- NA
hist(wht.stats$Fst)

# snps 
# set negative FSTs to NA
# this is sketchy, but for now i'm doing it
wht.stats.snp[is.infinite(wht.stats.snp$Fst) | wht.stats.snp$Fst<0, ]$Fst <- NA

wht.stats.snp %<>%
  filter(!is.na(Fst)) %>%
  mutate(fst.outlier=is.outlier(Fst,0.95,FALSE))

hist(wht.stats.snp$Fst)
max(wht.stats.snp$Fst)

mean(wht.stats.snp$Fst)


#### magic janky outlier detection with dplyr
is.outlier<-function(x, cutoff, equal.to){
  if (equal.to){
    return(x>=quantile(x,na.rm=TRUE,probs=cutoff)[1])
  }else{
    return(x>quantile(x,na.rm=TRUE,probs=cutoff)[1])
  }
}

wht.stats %<>%
  mutate(fst.outlier=is.outlier(Fst,0.99,FALSE))

wht.stats %>%
  filter(grepl("group",Chr))%>%
  ggplot(aes(x=StartPos, y=Dxy))+
  geom_point()+
  facet_grid(.~Chr)

wht.stats %>%
  filter(!is.na(Fst)) %>%
  ggplot(aes(x=Fst))+
  geom_histogram()+
  facet_grid(.~Chr)

wht.stats.snp %>%
  #filter(grepl("group",CHROM))%>%
  mutate(on.scaff = !grepl("group",CHROM)) %>%
  ggplot(aes(x=Fst, fill=fst.outlier))+
  #ggplot(aes(x=POS, y=Fst, color=fst.outlier))+
  geom_histogram()+
  facet_wrap(~on.scaff, scales = "free")
  