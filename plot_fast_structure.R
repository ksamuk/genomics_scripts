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

# assign membership

membership<-c()
for (i in 1:length(meanq.df$id)){
  k.value <- meanq.df[i,]$k.value.run
  member.columns <- meanq.df[i,4:(4+(k.value-1))]
  member.cluster<-as.numeric(which.max(member.columns))
  membership <- c(membership, member.cluster)
}

member.data <- data.frame(meanq.df, membership)
write.table(member.data, "whtcmn_gb2015_membership.txt")

meanq.df <- melt(meanq.df, 
          id.vars=c("pop","id","k.value.run"),
          variable.name = "k",
          value.name="q.value")

# reloaded to allow plyr/dplyr swap
library("dplyr")

# filter NAs
meanq.df <- meanq.df %>% 
  filter(!is.na(q.value))


# DISTRUCT-type plot


meanq.df %>%
  filter(k.value.run==3) %>%
  group_by(id) %>%
  

meanq.df%>%
  filter(k.value.run==2) %>%
  arrange(pop)%>%
  ggplot(aes(x=id, y=q.value, fill=factor(k)))+
  #ggplot(aes(x=id, y=q.value, fill=factor(pop)))+
    geom_bar(stat="identity", width=1)+
    #geom_bar(aes(x=id,fill=pop, y=.01),stat="identity",width=1,position="stack")+
    #geom_text(aes(label=pop))+
    theme_classic()+
    theme(axis.text.x=element_blank(), 
          axis.ticks=element_blank(), 
          axis.line=element_blank(),
          axis.title=element_blank())+
    facet_wrap(~pop, scales="free")

meanq.df%>%
  filter(k.value.run==3) %>%
  group_by(pop,id) %>%
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

meanq.df%>%
  filter(k.value.run==2) %>%
  group_by(pop,id) %>%
  mutate(q.order=max(q.value)
    
  