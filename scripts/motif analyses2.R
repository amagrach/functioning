

Now look at motifs and indirect interactions. To this end, we compare the roles of plant species at each site with low (<0.5) and high (>0.5) fruit set. We calculate for each plant species at each site the sum-normalised roles of all plant species in motifs up to 5 nodes and plot them using NMDS. ALthough we find a large overlap between the roles of both types ofplant species, we do find that plant species with high fruit set tend to be located in the higher values of NMDS1 which is related to network positions in which a generalist plant species interacts with highly specialized pollinator partners (e.g, positions 16 and 46, see Fig. 1 in bmotif paper, Simmons et al preprint)

#IB edited this in paper.Rmd
#First attempt to load the data, etc...

#prepare stuff-----
#If you don't have devtools install it
#install.packages("devtools")
library(devtools)
#Then install the data package. Reprot issues 
#install_github("ibartomeus/BeeFunData", ref = "fun")
library(BeeFunData)


# load packages
library(bmotif)
library(vegan)

#make the data available
#analysis...----


d<-all_interactions
f<-fruitset
s<-seedset


d2<-subset(d, d$Plant_gen_sp!="NA NA")
#d2$Plant_gen_sp

sites <- unique(d2$Site_ID)


empty_web <- function(x) {
  # Removes all-zero rows and columns from a matrix
  cempty <- which(colSums(x) == 0)
  rempty <- which(rowSums(x) == 0)
  if(length(rempty) != 0 & length(cempty) == 0){
    return(x[-rempty,,drop = FALSE])
  } else if(length(rempty) == 0 & length(cempty) != 0){
    return(x[,-cempty,drop = FALSE])
  } else if(length(rempty) != 0 & length(cempty) != 0){
    return(x[-rempty,-cempty,drop = FALSE])
  } else if(length(rempty) == 0 & length(cempty) == 0){
    return(x)
  }
}


webs <- list()
webs.or3<-list()
for(i in 1:length(sites)){
  temp <- subset(d2, Site_ID == sites[i])
  # for(j in 1:length(rounds)){
  #   temp2 <- subset(temp, Round == rounds[j])
  #   if (nrow(temp2) == 0) next
  
  web <- as.matrix.data.frame(table(temp$Plant_gen_sp, temp$Pollinator_gen_sp))
  web.or<-table(temp$Plant_gen_sp, temp$Pollinator_gen_sp)#keeps names of species
  
  web2<-empty_web(web)
  web.or2<-empty_web(web.or)
  
  webs.or3[[i]]<-web.or2 #keeps names of species
  webs[[i]]<-as.matrix(web2)
}

names(webs)<-sites
names(webs.or3)<-sites

out <- data.frame(web = NA, nrow = NA)

#calculates number of rows in each web and then sums everything to see how many rows we need in total
for(i in 1:length(webs)){
  s<- sum(nrow(webs[[i]]))
  
  n<-nrow(out)
  out[n + 1,1] <- as.character(sites[i])
  out[n + 1,2] <- s
  
}

out<-out[-1,]
out2<-sum(out$nrow)

mp <- as.data.frame(matrix(nrow = out2, ncol = 48, dimnames = list(NULL,c("web","species",paste0("m",1:46)))))
wr <- 0
for(i in names(webs)){
  
  pv <- positions(webs[[i]], six_node = FALSE, level = "rows", normalisation = "sum") # calculate species positions
  start_row <- which(is.na(mp$web))[1]
  end_row <- (start_row + nrow(pv)) - 1
  wr <- start_row:end_row
  mp[wr,"web"] <- i
  mp[wr,"species"] <- rownames(webs.or3[[i]])
  mp[wr,paste0("m",1:46)] <- pv
}



#SUBSET ONLY PLANTS FOR WHICH WE HAVE DATA ON FRUISET
#f4
f4$j<-paste(f4$Site_ID, f4$Plant.sp)

mp$j<-paste(mp$web, mp$species)


mp$mean.fruitset<-f4$mean.fruitset[match(mp$j, f4$j)]

mp2<-subset(mp, mp$mean.fruitset!="NA")
mp2$class[mp2$mean.fruitset>0.5]<-"high"
mp2$class[mp2$mean.fruitset<0.5]<-"low"



# How many native and introduced species?
mp2 <- mp2[order(mp2$class, decreasing = TRUE),]
class_f <- table(mp2$class)



# NMDS
example_NMDS <- metaMDS(vegdist(mp2[,paste0("m",1:46)], method = "bray"), k=2, distance = "bray", trymax = 100)
mds.fig <- ordiplot(example_NMDS, type = "none")
colors=c(rep("darkorange2",class_f["low"]),rep("purple",class_f["high"]))
treat <- mp2$class
for(i in unique(treat)) {
  ordihull(example_NMDS$point[grep(i,treat),],draw="polygon",
           groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } 
with(mp2, legend("topright", legend = unique(mp2$class), bty = "n",
                 col = unique(colors), pch = 21, pt.bg = unique(colors)))
points(mds.fig, "sites", pch = 19, col = "darkorange3", select = mp2$class == 
         "low")
points(mds.fig, "sites", pch = 19, col = "purple", select = mp2$class == 
         "high")

#example_NMDS <- metaMDS(mp2[,paste0("m",1:46)], k=2, distance = "bray", trymax = 100)
#example_NMDS$species # which positions are associated with which axes
#example_NMDS$species[order(example_NMDS$species[,"MDS2"]),]


#permanova

adonis(mp2[,3:48]~ mean.fruitset / web, data=mp2, method = "euclidean")

#permanova... 

adonis(mp2[,3:48]~ mean.fruitset, data=mp2, method = "euclidean", strata = mp2$web)
adonis(mp2[,3:48]~ class, data=mp2, method = "bray")




