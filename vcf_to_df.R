

source("http://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
library("VariantAnnotation")
library("dplyr")
library("magrittr")
library("ggplot2")
library("hierfstat")
library("data.table")

# read vcf file to 1/1 snp table

vcf.file <- list.files(pattern=".vcf")
whtcmn.gbs <- readGT(vcf.file)

# recode this useless format to data frame

lg <- rownames(whtcmn.gbs) %>% 
  strsplit(split=":") %>% 
  sapply(`[[`,1) %>% 
  gsub("group","",.) %>%
  as.roman %>%
  as.numeric

pos <- rownames(whtcmn.gbs) %>% 
  gsub("[^0-9]","",.) %>% 
  as.numeric

ref <- rownames(whtcmn.gbs) %>%
  strsplit(split="_") %>%
  sapply(`[[`,2) %>%
  strsplit(split="/") %>%
  sapply(`[[`,1)

alt <- rownames(whtcmn.gbs) %>%
  strsplit(split="_") %>%
  sapply(`[[`,2) %>%
  strsplit(split="/") %>%
  sapply(`[[`,2)

geno1 <- as.list(whtcmn.gbs) %>%
  gsub("[.]","NA/NA",.) %>%
  strsplit(split="/") %>% 
  sapply(`[`,1)

geno2 <- as.list(whtcmn.gbs) %>%
  gsub("[.]","NA/NA",.) %>%
  strsplit(split="/") %>% 
  sapply(`[`,2)

# 6 = number of underscores in id prefix
# "whtstbk_gbs_2015_L2_brds_SR80" to "SR80"
id <- colnames(whtcmn.gbs) %>%
  strsplit(split="_") %>%
  sapply(`[[`,6) %>%
  rep(each=length(whtcmn.gbs[,1]))

pop <- id %>%
  gsub("\\d","",.) %>%
  gsub("[:punct:]","",.)

whtcmn.df <- data.frame(pop, id, lg, pos, ref, alt, geno1, geno2)

whtcmn.df %<>% arrange(lg,pop,id,pos)

#fix factor NAs
whtcmn.df$geno1[whtcmn.df$geno1=="NA"] <- NA
whtcmn.df$geno2[whtcmn.df$geno2=="NA"] <- NA

whtcmn.df$geno1 %<>% as.character %>% as.numeric
whtcmn.df$geno2 %<>% as.character %>% as.numeric

whtcmn.df$ref.count <- as.numeric(whtcmn.df$geno1==0) + as.numeric(whtcmn.df$geno2==0)
whtcmn.df$alt.count <- as.numeric(whtcmn.df$geno1==1) + as.numeric(whtcmn.df$geno2==1)

# aaaah..that's better.

rm(list=c("whtcmn.gbs","geno1","geno2","id","pop","alt","ref","pos","lg"))
gc()

#assign genotypes
build.geno <- function (x){
  if(is.na(x[7])|is.na(x[8])){
    return("NN")
  }else{
    return(paste0(rep(x[5],x[9]),rep(x[6],x[10]),collapse=""))
  }
}

whtcmn.df$genotype <- rep("NN",length(whtcmn.df$genotype))
whtcmn.df$genotype <- apply(whtcmn.df[1:10,], MARGIN=1, build.geno)

whtcmn.df.fixed %>%
  


apply(whtcmn.df[1:10,], MARGIN=1, build.geno)

# write geno.df to file

gz1 <- gzfile("whtcmn.df.gz", "w")
write.csv(whtcmn.df, gz1, row.names = FALSE)
close(gz1)

whtcmn.df <- read.csv("whtcmn.df.gz", row.names = NULL)

#output "fake" FASTA file

df.to.fasta2 <- function (pop, id , genotype) {
  header <- paste0(">",pop[1],"_",id[1],"\n")
  dna <- paste0(genotype,collapse="")
  return(paste0(header,dna))
}

data.table(whtcmn.df)[, list(fasta=df.to.fasta2(pop,id,genotype)), by=key(whtcmn.df)]$fasta %>%
writeLines(file("test.fasta"))

whtcmn.df %>%
  filter(pop=="AL") %>%
  group_by(id)%>%
  summarise_each(funs(df.to.fasta)) %>%
  writeLines(file("test.fasta"))

# read in structure classes

structure.dat <- read.table("whtcmn_gb2015_membership.txt")
structure.dat %<>% filter(k.value.run==2) %>% select(pop,id,membership)

whtcmn.df <- left_join(whtcmn.df,structure.dat)

whtcmn.df %>%
  filter(lg==1)%>%
  filter(id=="AL1")%>%
  mutate(genotype=(geno1+geno2)) %>%
  ggplot(aes(x=pos,color=factor(genotype),y=1))+
  geom_bar(stat="identity")+
  scale_color_manual(values=c("blue","red","purple"))

whtcmn.df %>%
  filter(!is.na(geno1)|!is.na(geno2))%>%
  filter(!is.na(membership))%>%
  group_by(membership,lg,pos) %>%
  summarise(alt.frq=mean(alt.count)) %>% 
  ungroup %>%
  ungroup %>%
  ggplot(aes(x=pos,y=alt.frq,color=factor(membership)))+
  geom_smooth()+
  facet_wrap(~lg)

# hierstat format

whtcmn.wc <- t(whtcmn.gbs)
whtcmn.wc <- whtcmn.wc %>% 
  gsub("1/1","22",.) %>% 
  gsub("0/0","11",.) %>% 
  gsub("0/1","12",.) %>%
  gsub("[.]",NA,.) %>%
  as.data.frame

whtcmn.wc %<>% as.list %>% lapply(as.character) %>% lapply(as.numeric) %>% data.frame

# ids

id <- rownames(t(whtcmn.gbs)) %>%
  strsplit(split="_") %>%
  sapply(`[[`,6)

pop <- id %>%
  gsub("\\d","",.) %>%
  gsub("[:punct:]","",.)

# coerce and add cluster data
whtcmn.wc <- cbind(id,pop,whtcmn.wc)
whtcmn.wc <- left_join(whtcmn.wc,structure.dat)

#chop id and pop (implicit)
whtcmn.wc <-data.frame(whtcmn.wc [,length(whtcmn.wc)],whtcmn.wc [,-c(1,2,length(whtcmn.wc))])
names(whtcmn.wc)[1] <- "pop"

#calc wc

#had to manually do this part..
# write to file, go in and fix header
whtcmn.wc <- write.fstat("whtcmn.dat")
whtcmn.fst <- read.fstat("whtcmn.dat")

whtcmn.fst %<>%
  filter(!is.na(Pop))

whtcmn.fst.calc <- wc(whtcmn.fst)
#[1] "call"      "sigma"     "sigma.loc" "per.al"    "per.loc"   "FST"       "FIS"   

#overall FST: 0.015
whtcmn.fst.calc$FST

#back to df...
lg <- colnames(whtcmn.fst[,-1]) %>% 
  strsplit(split="[.]") %>% 
  sapply(`[[`,1) %>% 
  gsub("group","",.) %>%
  as.roman %>%
  as.numeric

pos <- colnames(whtcmn.fst[,-1]) %>% 
  gsub("[^0-9]","",.) %>% 
  as.numeric

fst <- whtcmn.fst.calc$per.loc$FST

whtcmn.fst.df <- data.frame(lg,pos,fst)
whtcmn.fst.df$fst[whtcmn.fst.df$fst<0]<-NA

# call outliers
is.outlier<-function(x){
  return(x>quantile(x,na.rm=TRUE,probs=0.95)[1])
}

whtcmn.fst.df <- whtcmn.fst.df%>%
  mutate(fst.outlier = is.outlier(fst))

# plots
whtcmn.fst.df %>%
  ggplot(aes(x=pos,y=fst, color=fst.outlier))+
  geom_point()+
  facet_wrap(~lg)

whtcmn.fst.df %>%
  ggplot(aes(x=pos,y=fst))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~lg)

whtcmn.fst.df %>%
  ggplot(aes(x=pos,y=as.numeric(fst.outlier)))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~lg)
