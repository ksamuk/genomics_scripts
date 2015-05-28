# build population map

install.packages(c("maps","mapproj","mapdata","rgeos","maptools","sp","raster","rgdal"))

library(maps)
library(mapproj)
library(mapdata)
library(rgeos)
library(maptools)
library(sp)
library(raster)
library(rgdal)

can1 <- getData('GADM', country="CAN", level=1)

# plot nova soctia            
#-66.39584, -59.6425, 43.38889, 47.22745  (xmin, xmax, ymin, ymax)

#zoomed coords
x.zoom <- c(-66.39584, -59.6425)
y.zoom <- c(43.38889, 47.22745)
  
plot(can1[7,], xlim=x.zoom, ylim= y.zoom)

# Extent rectangle for inset map


p1<-ggplot()+
  geom_polygon(data=can1[7,], aes(long,lat, group=group), color="black",fill="grey50")+
  coord_map(xlim = x.zoom ,ylim = y.zoom)+
  theme_classic()+
  theme(axis.line=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank())

# Bras d'Or extents
x.zoom2 <- c(-61.20, -60.20)
y.zoom2 <- c(45.60, 46.35)
plot(can1[7,], xlim=x.zoom2, ylim= y.zoom2)
    
p2<-ggplot()+
  #geom_rect(aes(xmax=x.zoom2[2],xmin=x.zoom2[1],ymax=y.zoom2[2],ymin=y.zoom2[1], fill="skyblue"))+
  geom_polygon(data=can1[7,], aes(long,lat, group=group), color="black",fill="grey50")+
  coord_map(xlim = x.zoom2 ,ylim = y.zoom2)+
  theme_classic()+
  theme(axis.line=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank())
  
  png(file="mrdq.png",w=1800,h=1800, res=300)
  grid.newpage()
  v2<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
  v1<-viewport(width = 0.3, height = 0.3, x = 0.86, y = 0.28) #plot area for the inset map
  print(p1,vp=v2) 
  print(p2,vp=v2)
  dev.off()
