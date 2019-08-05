
library(BeeFunData)

si<-BeeFunData::sites #lat long data for sites 

f.agg$lat<-si$latitude[match(f.agg$Site_ID, si$Site_ID)]
f.agg$lon<-si$longitude[match(f.agg$Site_ID, si$Site_ID)]

s.comm$lat<-si$latitude[match(s.comm$Site_ID, si$Site_ID)]
s.comm$lon<-si$longitude[match(s.comm$Site_ID, si$Site_ID)]

library(ecodist)

#distance matrices
xy <- cbind(as.vector(f.agg$lat), as.vector(f.agg$long))
xydis <-distance(xy, method="euclidean")

nesteddis<-distance(as.vector(f.agg$nodf.song), method="euclidean")
functdis<-distance(as.vector(f.agg$functional.comp.poll), method="euclidean")
richdis<-distance(as.vector(f.agg$species.poll), method="euclidean")
fruitdis<-distance(as.vector(f.agg$mean.fruitset), method="euclidean")
visdis<-distance(as.vector(f.agg$tot.visits), method="euclidean")
seeddis<-distance(as.vector(s.comm$mean.seedset), method="euclidean")
seedwdis<-distance(as.vector(s.comm$mean.seedweight), method="euclidean")
fruitwdis<-distance(as.vector(s.comm$mean.fruitweight), method="euclidean")

#relationship among matrices
Nested_dist <- mantel(xydis ~ nesteddis, nperm=1000)
Funct_dist <- mantel(xydis ~ functdis, nperm=1000)
Rich_dist <- mantel(richdis ~ xydis, nperm=1000)
Vis_dist <- mantel(visdis ~ xydis, nperm=1000)
Fruit_dist <- mantel(fruitdis ~ xydis, nperm=1000)
Seed_dist <- mantel(seeddis ~ xydis, nperm=1000)
Seedw_dist <- mantel(seedwdis ~ xydis, nperm=1000)
Fruitw_dist <- mantel(fruitwdis ~ xydis, nperm=1000)

Varmantel <- cbind(Nested_dist,Funct_dist, Rich_dist, Vis_dist, Fruit_dist, Seed_dist, Seedw_dist, Fruitw_dist)
Varmantel


#correlograms
Nestcor=mgram(nesteddis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Nestcor, main="Nestedness")

Functcor=mgram(functdis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Functcor, main="Pollinator niche complementarity")

Richcor=mgram(richdis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Richcor, main="Pollinator species richness")

Viscor=mgram(visdis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Viscor, main="Total visits")

Fruitcor=mgram(fruitdis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Fruitcor, main="Fruitset")

Seedcor=mgram(seeddis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Seedcor, main="Seeds per fruit")

Seedwcor=mgram(seedwdis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Seedwcor, main="Seed weight")

Fruitwcor=mgram(fruitwdis, xydis, nperm=10000, nboot=100, nclass=7)
plot(Fruitwcor, main="Fruit weight")

par(mfrow=c(3,3))
plot(Nestcor, main="Nestedness");plot(Functcor, main="Pollinator niche complementarity");
plot(Richcor, main="Pollinator species richness");plot(Viscor, main="Total visits");
plot(Fruitcor, main="Fruit set"); plot(Seedcor, main="Seeds per fruit"); plot(Seedwcor, main="Seed weight");
plot(Fruitwcor, main="Fruit weight")




