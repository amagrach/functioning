
#library(devtools)
#Then install the data package. Reprot issues 
#install_github("ibartomeus/BeeFunData", ref = "fun")

library(BeeFunData)
require(bipartite)
source('toolbox.R') #load the toolbox from Song et al paper
#devtools::install_github("CHoeppke/maxnodf")
library(maxnodf)#load improved algorithm from Simmons et al
library(plyr)

d<-all_interactions



#Repeat analyses removing out of transect observations of rare specimens

d.sinout <- subset(d, !Out %in% c("out"))



sites.sinout <- unique(d.sinout$Site_ID)


out.site.sinout <- data.frame(Site_id = NA,  
                              species.poll=NA, species.pl=NA, 
                              functional.comp.poll=NA, functional.comp.plant=NA,
                              nodf.song=NA)

outsp.site.sinout <- data.frame(Site_id = NA, species=NA,  norm_degree = NA,
                                weigh_closeness=NA, d=NA)


outsp.pl.site.sinout <- data.frame(Site_id = NA, species=NA,  norm_degree = NA,
                                   weigh_closeness=NA, d=NA)


webs <- list()
for(i in 1:length(sites.sinout)){
  print(i)
  temp <- subset(d.sinout, Site_ID == sites.sinout[i])
  
  
  web<-table(temp$Plant_gen_sp, temp$Pollinator_gen_sp)
  
  #p1<-plotweb(web)
  spntw <- try(specieslevel(web), TRUE)
  ntw <- networklevel(web)
  if(isTRUE(class(spntw)=="try-error")) {next} 
  
  #song et al nestedness corrected with NODF max by Simmons et al
  
  #remove all columns and rows with all values=0
  web3 <- web[, colSums(web != 0) > 0] 
  web4 <- web3[rowSums(web3[, -1] > 0) != 0, ]
  
  NODF <- nestedness_NODF(web4) # this calculates the raw value of NODF as in Song
  if(i %in% c(6,15,16)){ #For sites wchich don't comply, use Song. Only three sites.
    max_NODF <- max_nest(web4[,-1])
    combined_NODF <- comb_nest(web4[,-1],NODF,max_NODF) # this calculates the combined NODF statistic as described in the manuscript
  } else {
    max_NODF <- maxnodf(web=web4, quality=2) # this calculates the maximum value of NODF for that network, based on SImmons correction of Song
    print(paste("Simmons = ", max_NODF$max_nodf, ", Song = ", max_nest(web4)))
    combined_NODF <- comb_nest(web4,NODF,max_NODF$max_nodf) # this calculates the combined NODF statistic as described in the manuscript
  }
  
  
  n <- nrow(out.site.sinout)
  webs[[n]] <- web4  
  n2 <- nrow(outsp.site.sinout)
  n3 <- nrow(outsp.pl.site.sinout)
  
  out.site.sinout[n + 1,1] <- as.character(sites[i])
  out.site.sinout[n + 1,2:5] <- c(ntw[20:21], ntw[42:43])
  out.site.sinout[n + 1,6] <- try(combined_NODF, TRUE)
  
  outsp.site.sinout[n2+seq(nrow(spntw$`higher level`)),1] <- as.character(sites[i])
  outsp.site.sinout[n2+seq(nrow(spntw$`higher level`)),2] <- try(as.character(rownames(spntw$`higher level`)), TRUE)
  outsp.site.sinout[n2+seq(nrow(spntw$`higher level`)),3:5] <- try(c(spntw$`higher level`[2], 
                                                                     spntw$`higher level`[14],
                                                                     spntw$`higher level`[20]), TRUE)
  
  
  
  outsp.pl.site.sinout[n3+seq(nrow(spntw$`lower level`)),1] <- as.character(sites[i])
  outsp.pl.site.sinout[n3+seq(nrow(spntw$`lower level`)),2] <- try(as.character(rownames(spntw$`lower level`)), TRUE)
  outsp.pl.site.sinout[n3+seq(nrow(spntw$`lower level`)),3:5] <- try(c(spntw$`lower level`[2], 
                                                                       spntw$`lower level`[14],
                                                                       spntw$`lower level`[20]), TRUE)
  
}  
str(out.site.sinout)
head(out.site.sinout)
out.site.sinout

head(outsp.site.sinout)
str(outsp.site.sinout)

head(outsp.pl.site.sinout)
str(outsp.pl.site.sinout)

write.table(out.site.sinout[2:nrow(out.site.sinout),], file= "SITE_network_level_metrics_sinout.csv", row.names= FALSE, sep= ",")
write.table(outsp.site.sinout[2:nrow(outsp.site.sinout),], file= "SITE_species_level_metrics_sinout.csv", row.names= FALSE, sep= ",")
write.table(outsp.pl.site.sinout[2:nrow(outsp.pl.site.sinout),], file= "SITE_plant_species_level_metrics_sinout.csv", row.names= FALSE, sep= ",")



##network metrics removing out of transect observations
outsp.pl.site.sinout<-read.csv("SITE_plant_species_level_metrics_sinout.csv")
outsp.site.sinout<-read.csv("SITE_species_level_metrics_sinout.csv")
out.site.sinout<-read.csv("SITE_network_level_metrics_sinout.csv")


# CALCULATE NICHE OVERLAP BETWEEN PLANT SPECIES AND TOTAL NUMBER OF VISITS AND JOIN FRUITSET DATA WITH NETWORK METRICS

#calculate total number of visits per plant sps

d.vis.sinout<-ddply(d.sinout, c("Site_ID", "Plant_gen_sp"), summarise,
             tot.visits= sum(Frequency))



d.vis.sinout$j<-paste(d.vis.sinout$Site_ID, d.vis.sinout$Plant_gen_sp)

#calculate average niche overlap per species 



sites <- unique(d.vis.sinout$Site_ID)

out.niche.sinout <- data.frame(Site_id = NA, row = NA, col = NA,value  = NA)

webs <- list()
for(i in 1:length(sites)){
  temp <- subset(d.sinout, Site_ID == sites[i])
  temp2<-droplevels(temp)
  
  
  web<-table(temp2$Pollinator_gen_sp, temp2$Plant_gen_sp)
  

  ni<-niche.overlap(web, method="morisita") #THIS SHOULD BE EXPLAINED IN METHODS
  df <- melt(as.matrix(ni), varnames = c("row", "col"))
  df$row<-as.character(df$row)
  df$col<-as.character(df$col)
  
  n <- nrow(out.niche.sinout)
  
  out.niche.sinout[n + seq(nrow(df)),1] <- as.character(sites[i])
  out.niche.sinout[n + seq(nrow(df)),2:3] <- (df[,1:2])
  out.niche.sinout[n + seq(nrow(df)),4] <-df[,3]
  
}

out.niche.agg.sinout <- ddply(out.niche.sinout, c("Site_id", "col"), summarise,
                       mean.niche.over    = mean(value))


out.niche.agg.sinout$j<-paste(out.niche.agg.sinout$Site_id, out.niche.agg.sinout$col)

#use unaggregated data and join it with network metrics

f2$functional.comp.poll.sinout<-out.site.sinout$functional.comp.poll[match(f2$Site_ID, out.site.sinout$Site_id)]
f2$functional.comp.pl.sinout<-out.site.sinout$functional.comp.pl[match(f2$Site_ID, out.site.sinout$Site_id)]
f2$nodf.song.sinout<-out.site.sinout$nodf.song[match(f2$Site_ID, out.site.sinout$Site_id)]
f2$species.poll.sinout<-out.site.sinout$species.poll[match(f2$Site_ID, out.site.sinout$Site_id)]

#also match with plant species level network data
f2$match<-paste(f2$Site_ID, f2$Plant.sp)
outsp.pl.site.sinout$match<-paste(outsp.pl.site.sinout$Site_id, outsp.pl.site.sinout$species)

f2$norm_degree.sinout<-outsp.pl.site.sinout$norm_degree[match(f2$match, outsp.pl.site.sinout$match)]
f2$weigh_closeness.sinout<-outsp.pl.site.sinout$weigh_closeness[match(f2$match, outsp.pl.site.sinout$match)]

#and with niche overlap and total number of visits
f2$niche.overlap.sinout<-out.niche.agg.sinout$mean.niche.over[match(f2$match, out.niche.agg.sinout$j)]

f2$tot.visits.sinout<-d.vis.sinout$tot.visits[match(f2$match, d.vis.sinout$j)]
f2$tot.visits.sinout[is.na(f2$tot.visits.sinout)]<-0


#calculate number of pollinator species per plant species instead of norm_degree

sites <- unique(d.sinout$Site_ID)


out.sinout <- data.frame(Site_id = NA, plant.sps = NA, poll.sps = NA)


for(i in 1:length(sites)){
  temp <- subset(d.sinout, Site_ID == sites[i])
  
  plants<-unique(temp$Plant_gen_sp)
  for(j in 1:length(plants)){
    temp2 <- subset(temp, Plant_gen_sp == plants[j])
    temp3<-droplevels(temp2)
    
    poll<-nlevels(unique(temp3$Pollinator_gen_sp))
    
    n <- nrow(out.sinout)
    
    out.sinout[n + 1,1] <- as.character(sites[i])
    out.sinout[n + 1,2] <- as.character(plants[j])
    out.sinout[n + 1,3] <- poll
    
  }
}
#out.sinout


out.sinout2<-out.sinout[-1,]
out.sinout2$match<-paste(out.sinout2$Site_id, out.sinout2$plant.sps)

f2$pollspsperplant.sinout<-out.sinout2$poll.sps[match(f2$match, out.sinout2$match)]


#scale all variables used

f2$functional.comp.poll.sinout.sc<-scale(f2$functional.comp.poll.sinout, center = T, scale = T)
f2$functional.comp.pl.sinout.sc<-scale(f2$functional.comp.pl.sinout, center = T, scale = T)
f2$nodf.song.sinout.sc<-scale(f2$nodf.song.sinout, center = T, scale = T)
f2$species.poll.sinout.sc<-scale(f2$species.poll.sinout, center = T, scale = T)
f2$norm_degree.sinout.sc<-scale(f2$norm_degree.sinout, center = T, scale = T)
f2$weigh_closeness.sinout.sc<-scale(f2$weigh_closeness.sinout, center = T, scale = T)
f2$niche.overlap.sinout.sc<-scale(f2$niche.overlap.sinout, center = T, scale = T)
f2$pollspsperplant.sinout.sc<-scale(f2$pollspsperplant.sinout, center = T, scale = T)
f2$tot.visits.sinout.sc<-scale(f2$tot.visits.sinout, center = T, scale = T)



#compare 2 models: one with simple diversity metrics (pollinator richness and number of visits), compared to one where position in a network (centrality) and niche overlap are included with these variables.

m1.sinout<-glmer(tot.fruitset ~  pollspsperplant.sinout.sc + tot.visits.sinout.sc  + (1|Plant.sp:Site_ID) + (1|Site_ID), 
                 family="binomial", data=f2)

summary(m1.sinout)
vif(m1.sinout)  #check for correlations between variables included in model



m2.sinout<-glmer(tot.fruitset ~  pollspsperplant.sinout.sc  + tot.visits.sinout.sc + weigh_closeness.sinout.sc + 
                   niche.overlap.sinout.sc + (1|Plant.sp:Site_ID)+ (1|Site_ID), family="binomial", data=f2)

summary(m2.sinout) 
vif(m2.sinout)

#select best model based on AIC 

AICtab(m1.sinout, m2.sinout) 

r2<-rsquared(c(m2.sinout, m1.sinout)) #check R^2 values for both models 

# extract estimates and st errors to present in table
coefs.s11A <- data.frame(coef(summary(m2.sinout)))
TableS11A<-round(coefs.s11A[,1:3],digits=2)
rownames(TableS11A)<-c("(Intercept)", "Pollinator richness", "Total number of visits", "Centrality", "Plant niche overlap")



##seed set

s3$functional.comp.poll.sinout<-out.site.sinout$functional.comp.poll[match(s3$Site_ID, out.site.sinout$Site_id)]
s3$functional.comp.pl.sinout<-out.site.sinout$functional.comp.pl[match(s3$Site_ID, out.site.sinout$Site_id)]
s3$nodf.song.sinout<-out.site.sinout$nodf.song[match(s3$Site_ID, out.site.sinout$Site_id)]
s3$species.poll.sinout<-out.site.sinout$species.poll[match(s3$Site_ID, out.site.sinout$Site_id)]

#also match with plant species level network data
s3$match<-paste(s3$Site_ID, s3$Plant.sp)
outsp.pl.site.sinout$match<-paste(outsp.pl.site.sinout$Site_id, outsp.pl.site.sinout$species)

s3$norm_degree.sinout<-outsp.pl.site.sinout$norm_degree[match(s3$match, outsp.pl.site.sinout$match)]
s3$weigh_closeness.sinout<-outsp.pl.site.sinout$weigh_closeness[match(s3$match, outsp.pl.site.sinout$match)]
s3$niche.overlap.sinout<-out.niche.agg.sinout$mean.niche.over[match(s3$match, out.niche.agg.sinout$j)]


s3$tot.visits.sinout<-d.vis.sinout$tot.visits[match(s3$match, d.vis.sinout$j)]
s3$tot.visits.sinout[is.na(s3$tot.visits.sinout)]<-0

s3$pollspsperplant.sinout<-out.sinout2$poll.sps[match(s3$match, out.sinout2$match)]


#scale all variables used

s3$norm_degree.sinout.sc<-scale(s3$norm_degree.sinout, center = T, scale = T)
s3$weigh_closeness.sinout.sc<-scale(s3$weigh_closeness.sinout, center = T, scale = T)
s3$niche.overlap.sinout.sc<-scale(s3$niche.overlap.sinout, center = T, scale = T)
s3$pollspsperplant.sinout.sc<-scale(s3$pollspsperplant.sinout, center = T, scale = T)
s3$tot.visits.sinout.sc<-scale(s3$tot.visits.sinout, center = T, scale = T)

s3$mean.seedn.sc<-scale(s3$mean.seedn, center = T, scale = T)


#fit two models

m1.seed.sinout<-glmer(mean.seedn.sc ~  pollspsperplant.sinout.sc + tot.visits.sinout.sc  + (1|Plant.sp:Site_ID) + (1|Site_ID), family="gaussian", data=s3)

summary(m1.seed.sinout)
vif(m1.seed.sinout)


m2.seed.sinout<-glmer(mean.seedn.sc ~  pollspsperplant.sinout.sc  + tot.visits.sinout.sc + 
                        weigh_closeness.sinout.sc + niche.overlap.sinout.sc  + (1|Plant.sp:Site_ID) + (1|Site_ID), 
                      family="gaussian", data=s3)


summary(m2.seed.sinout)
vif(m2.seed.sinout)


AICtab(m1.seed.sinout, m2.seed.sinout) #m2 is best

r2b<-rsquared(c(m2.seed, m1.seed)) #check R^2 values for both models 


# extract estimates and st errors to present in table
coefs.s11B <- data.frame(coef(summary(m2.seed.sinout)))
TableS11B<-round(coefs.s11B[,1:3],digits=2)
rownames(TableS11B)<-c("(Intercept)", "Pollinator richness", "Total number of visits", "Centrality", "Plant niche overlap")


#SEED WEIGHT

m1.seed.w.sinout<-glmmTMB(mean.seedw ~  pollspsperplant.sinout.sc + tot.visits.sinout.sc  + 
                            (1|Plant.sp:Site_ID) + (1|Site_ID), family="poisson", data=s3)

summary(m1.seed.w.sinout)


m2.seed.w.sinout<-glmmTMB(mean.seedw ~  pollspsperplant.sinout.sc + tot.visits.sinout.sc + weigh_closeness.sinout.sc + 
                            niche.overlap.sinout.sc  + 
                     (1|Plant.sp:Site_ID) + (1|Site_ID), family="poisson", data=s3)

summary(m2.seed.w.sinout)


#select best model based on AIC 

AICtab(m1.seed.w.sinout, m2.seed.w.sinout) 


# extract estimates and st error for m1 and present them in paper
tableS11C<-round(summary(m1.seed.w.sinout)$coefficients$cond[,1:3], digits=2)
rownames(tableS11C)<-c("(Intercept)", "Pollinator richness", "Total number of visits")


#FRUIT WEIGHT

m1.fruit.w.sinout<-glmmTMB(mean.fruitweight ~  pollspsperplant.sinout.sc + tot.visits.sinout.sc  + (1|Plant.sp:Site_ID) + (1|Site_ID), family="poisson", data=s3)

summary(m1.fruit.w.sinout)



m2.fruit.w.sinout<-glmmTMB(mean.fruitweight ~  pollspsperplant.sinout.sc  + tot.visits.sinout.sc + weigh_closeness.sinout.sc + 
                             niche.overlap.sinout.sc  + 
                      (1|Plant.sp:Site_ID) + (1|Site_ID), family="poisson", data=s3)

summary(m2.fruit.w.sinout)


#select best model based on AIC 

AICtab(m1.fruit.w.sinout, m2.fruit.w.sinout) 

# extract estimates and st error for m1 and present them in paper
tableS11D<-round(summary(m1.fruit.w.sinout)$coefficients$cond[,1:3], digits=2)
rownames(tableS11D)<-c("(Intercept)", "Pollinator richness", "Total number of visits")


######################## NTW STRUCTURE

#aggregate values at site level

f.ag.sinout <- ddply(f2, c("Site_ID", "Plant.sp", "nodf.song.sinout", "functional.comp.poll.sinout",
                    "functional.comp.pl.sinout", "species.poll.sinout"), 
              summarise, tot.visits.sinout    = sum(tot.visits.sinout),
              mean.fruitset = mean(tot.fruitset))

f.agg.sinout <- ddply(f.ag.sinout, c("Site_ID",  "nodf.song.sinout", "functional.comp.poll.sinout", 
                       "functional.comp.pl.sinout","species.poll.sinout"), 
               summarise, tot.visits.sinout    = sum(tot.visits.sinout),
               mean.fruitset = mean(mean.fruitset))


#scale all variables to be able to compare effect sizes 

f.agg.sinout$nodf.song.sinout.sc <- scale(f.agg.sinout$nodf.song.sinout, center = T, scale = T)
f.agg.sinout$functional.comp.poll.sinout.sc <- scale(f.agg.sinout$functional.comp.poll.sinout, center = T, scale = T)
f.agg.sinout$functional.comp.pl.sinout.sc <- scale(f.agg.sinout$functional.comp.pl.sinout, center = T, scale = T)
f.agg.sinout$species.poll.sinout.sc <- scale(f.agg.sinout$species.poll.sinout, center = T, scale = T)
f.agg.sinout$tot.visits.sinout.sc <- scale(f.agg.sinout$tot.visits.sinout, center = T, scale = T)

#FRUITset

m1.fruit.net.sinout<-glmmTMB(mean.fruitset ~  species.poll.sinout.sc + tot.visits.sinout.sc, family="beta_family", data=f.agg.sinout)

summary(m1.fruit.net.sinout)

m2.fruit.net.sinout<-glmmTMB(mean.fruitset ~  species.poll.sinout.sc + tot.visits.sinout.sc + nodf.song.sinout.sc + 
                               functional.comp.poll.sinout.sc, family="beta_family", data=f.agg.sinout)

summary(m2.fruit.net.sinout)


#select best model based on AIC 

AICtab(m1.fruit.net.sinout, m2.fruit.net.sinout)


#seed number


s.com.sinout <- ddply(s3, c("Site_ID", "Plant.sp", "nodf.song.sinout", "functional.comp.poll.sinout",
                     "functional.comp.pl.sinout", "species.poll.sinout"), 
               summarise, tot.visits.sinout    = sum(tot.visits.sinout),
               mean.seedset = mean(mean.seedn),
               mean.seedweight = mean(mean.seedw),
               mean.fruitweight = mean(mean.fruitweight))


s.comm.sinout <- ddply(s.com.sinout, c("Site_ID", "nodf.song.sinout", "functional.comp.poll.sinout", 
                         "functional.comp.pl.sinout","species.poll.sinout"), 
                summarise, tot.visits.sinout    = sum(tot.visits.sinout),
                mean.seedset = mean(mean.seedset),
                mean.seedweight = mean(mean.seedweight),
                mean.fruitweight = mean(mean.fruitweight))



#scale variables for comparison of effect sizes

s.comm.sinout$nodf.song.sinout.sc <- scale(s.comm.sinout$nodf.song.sinout, center = T, scale = T)
s.comm.sinout$functional.comp.poll.sinout.sc <- scale(s.comm.sinout$functional.comp.poll.sinout, center = T, scale = T)
s.comm.sinout$functional.comp.pl.sinout.sc <- scale(s.comm.sinout$functional.comp.pl.sinout, center = T, scale = T)
s.comm.sinout$species.poll.sinout.sc <- scale(s.comm.sinout$species.poll.sinout, center = T, scale = T)
s.comm.sinout$tot.visits.sinout.sc <- scale(s.comm.sinout$tot.visits.sinout, center = T, scale = T)

#NUMBER OF SEEDS PER FRUIT


m1.seed.net.sinout <- lm(mean.seedset ~  species.poll.sinout.sc + tot.visits.sinout.sc,
                  data=s.comm.sinout,family="gaussian")

summary(m1.seed.net.sinout)

m2.seed.net.sinout <- lm(mean.seedset ~  species.poll.sinout.sc + tot.visits.sinout.sc + nodf.song.sinout.sc + 
                           functional.comp.poll.sinout.sc, data=s.comm.sinout,family="gaussian")

summary(m2.seed.net.sinout)

#select best model based on AIC 
AICtab(m1.seed.net.sinout, m2.seed.net.sinout)

#FRUIT WEIGHT


m1.fruitw.net.sinout<-lm(mean.fruitweight ~  species.poll.sinout.sc + tot.visits.sinout.sc, family="gaussian", data=s.comm.sinout)

summary(m1.fruitw.net.sinout)

vif(m1.fruitw.net.sinout)

m2.fruitw.net.sinout<-lm(mean.fruitweight ~  species.poll.sinout.sc + tot.visits.sinout.sc + nodf.song.sinout.sc + 
                           functional.comp.poll.sinout.sc, family="gaussian", data=s.comm.sinout)

summary(m2.fruitw.net.sinout)

vif(m2.fruitw.net.sinout)

#select best model based on AIC 

AICtab(m1.fruitw.net.sinout, m2.fruitw.net.sinout)

#SEED WEIGHT


m1.seedw.net.sinout<-glm(mean.seedweight ~  species.poll.sinout.sc + tot.visits.sinout.sc, family="gaussian", data=s.comm.sinout)

summary(m1.seedw.net.sinout)


vif(m1.seedw.net.sinout)

m2.seedw.net.sinout<-glm(mean.seedweight ~  species.poll.sinout.sc + tot.visits.sinout.sc + nodf.song.sinout.sc + 
                           functional.comp.poll.sinout.sc, family="gaussian", data=s.comm.sinout)

summary(m2.seedw.net.sinout)

vif(m2.seedw.net.sinout)

#select best model based on AIC 

AICtab(m1.seedw.net.sinout, m2.seedw.net.sinout) #m1 is better

