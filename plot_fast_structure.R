#### plot fastStructure output with ggplot

library("ggplot2")
library("dplyr")
library("plyr")
library("reshape2")

# a directory containing sample id files (single column of IDs from ped file) and meanq files
output.dir <- "fastStructure"

# the file names
sample.id.file <- list.files(output.dir, pattern=".txt", full.names=TRUE)
meanq.files <- list.files(output.dir, pattern=".meanQ", full.names=TRUE)

# strip id prefixes
# assumes snake_case_ids e.g. gbs_2015_POP##

ids.raw <- scan(file=sample.id.file, what=character())

strip_prefix <- function(x){
  split <- strsplit(x, split="_")
  return(split[[1]][length(split[[1]])])
}

ids <- sapply(ids.raw, strip_prefix)
ids <- unname(ids)

# get pop ids from ids
# \\d = digits, \\W = any weird punct
pop.ids <- gsub("\\d", "", ids)
pop.ids <- gsub("\\W", "", pop.ids)

# build id/pop frame  

id.df <- data.frame(pop=pop.ids, id=ids)

# merge in q values

meanq.df <- data.frame()

for (i in 1:length(meanq.files)){
  
  # read meanq file
  q.file <- read.table(meanq.files[i], header=FALSE)
  
  # get k value
  k.value.run <- length(q.file)
  
  # fill in column names
  names(q.file) <- (1:k.value.run)
  
  # fill in k.value as column
  k.value.run <- rep(k.value.run, length(q.file[,1]))
  
  # bind to df
  q.file <- cbind(k.value.run, q.file)

  # bind to id.df
  q.file <- cbind(id.df, q.file)
  
  # bind with previous (if not first run )
  if (i==1){
    meanq.df <- q.file
  } else{
    meanq.df <- rbind.fill(meanq.df, q.file)
  }
  
}

meanq.df <- melt(meanq.df, 
          id.vars=c("pop","id","k.value.run"),
          variable.name = "k",
          value.name="q.value")

# reloaded to allow plyr/dplyr swap
library("dplyr")

# filter NAs
meanq.df <- meanq.df %>% 
  filter(!is.na(q.value))

meanq.df.exag <- meanq.df

meanq.df.exag[meanq.df.exag$k==3 & meanq.df.exag$k.value.run==3, ]$q.value <- meanq.df.exag[meanq.df.exag$k==3 & meanq.df.exag$k.value.run==3, ]$q.value

# DISTRUCT-type plot

meanq.df%>%
  filter(k.value.run==2) %>%
  arrange(pop)%>%
  ggplot(aes(x=id, y=q.value, fill=factor(k)))+
  #ggplot(aes(x=id, y=q.value, fill=factor(pop)))+
    geom_bar(stat="identity", width=1)+
    #geom_text(aes(label=pop))+
    theme_classic()+
    theme(axis.text.x=element_blank(), 
          axis.ticks=element_blank(), 
          axis.line=element_blank(),
          axis.title=element_blank())

meanq.df%>%
  filter(k.value.run==2) %>%
  group_by(pop) %>%
  mutate(q.max=max(q.value))%>%
  ungroup() %>%
  arrange(q.max)%>%
  ggplot(aes(x=id, y=q.value, fill=factor(k)))+
  #ggplot(aes(x=id, y=q.value, fill=factor(pop)))+
  geom_bar(stat="identity", width=1)+
  #geom_text(aes(label=pop))+
  theme_classic()+
  theme(axis.text.x=element_blank(), 
        axis.ticks=element_blank(), 
        axis.line=element_blank(),
        axis.title=element_blank())

assignment<-c()
for (i in 1:length(meanq.df$qvalue)){
  if meanq.df$qvalue >0.0
}
    
  