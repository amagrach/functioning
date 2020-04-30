
#library(devtools)
#Then install the data package. Reprot issues 
#install_github("ibartomeus/BeeFunData", ref = "fun")

library(BeeFunData)
require(bipartite)
source('toolbox.R') #load the toolbox from Song et al paper
#devtools::install_github("CHoeppke/maxnodf")
library(maxnodf)#load improved algorithm from Simmons et al
library(plyr)
library(nnet)

d<-all_interactions

dtrans<-subset(d, d$Out=="transect")
dtrans<-droplevels(dtrans)



##1st create a bipartite matrix for each site MERGING ALL ROUNDS TOGETHER and extract a set of metrics at the network level

sites <- unique(dtrans$Site_ID)


out.site <- data.frame(Site_id = NA,  
                       species.poll=NA, species.pl=NA, 
                       functional.comp.poll=NA, functional.comp.plant=NA, 
                       nodf.song=NA, nodf.song2=NA)

outsp.site <- data.frame(Site_id = NA, species=NA,  norm_degree = NA, nested.rank=NA, strength =NA,
                         weigh_closeness=NA, d=NA, nested_contribution=NA)


outsp.pl.site <- data.frame(Site_id = NA, species=NA,  norm_degree = NA, nested.rank=NA,strength =NA,
                            weigh_closeness=NA, d=NA, nested_contribution=NA)



webs <- list()
for(i in 1:length(sites)){
  print(i)
  temp <- subset(dtrans, Site_ID == sites[i])
  
  

  web<-table(temp$Plant_gen_sp, temp$Pollinator_gen_sp)
  
   #p1<-plotweb(web)
  spntw <- try(specieslevel(web), TRUE)
 
  ntw <- networklevel(web)
  if(isTRUE(class(spntw)=="try-error")) {next} 
  
  #song et al nestedness corrected with NODF max by Simmons et al
  
  #remove all columns and rows with all values=0
  web3 <- web[, colSums(web != 0) > 0] 
  web4 <- web3[rowSums(web3 > 0) != 0, ]
  
  nestedcontr<-nestedcontribution(web4, nsimul=99)
  
   NODF <- nestedness_NODF(web4) # this calculates the raw value of NODF as in Song
    if(i %in% c(1,2,4, 5,6,7,8, 10,11,12,14, 15,16)){ #For sites wchich don't comply, use Song. 
      max_NODF <- max_nest(web4[,-1])
      combined_NODF <- comb_nest(web4[,-1],NODF,max_NODF) # this calculates the combined NODF statistic as described in the manuscript
    } else {
      max_NODF <- maxnodf(web=web4, quality=2) # this calculates the maximum value of NODF for that network, based on SImmons correction of Song
      print(paste("Simmons = ", max_NODF$max_nodf, ", Song = ", max_nest(web4)))
      combined_NODF <- comb_nest(web4,NODF,max_NODF$max_nodf) # this calculates the combined NODF statistic as described in the manuscript
    }
  
   NODF2 <- nestedness_NODF(web4) # recalculate everything using only Song approach to see correlation between both approaches
   
   max_NODF2 <- max_nest(web4[,-1])
   combined_NODF2 <- comb_nest(web4[,-1],NODF2,max_NODF2)  
  
  n <- nrow(out.site)
  webs[[n]] <- web4  
  n2 <- nrow(outsp.site)
  n3 <- nrow(outsp.pl.site)
  
  out.site[n + 1,1] <- as.character(sites[i])
  out.site[n + 1,2:5] <- c(ntw[20:21], ntw[42:43])
  out.site[n + 1,6] <- try(combined_NODF, TRUE)
  out.site[n + 1,7] <- try(combined_NODF2, TRUE)
 
  outsp.site[n2+seq(nrow(spntw$`higher level`)),1] <- as.character(sites[i])
  outsp.site[n2+seq(nrow(spntw$`higher level`)),2] <- try(as.character(rownames(spntw$`higher level`)), TRUE)
  outsp.site[n2+seq(nrow(spntw$`higher level`)),3:7] <- try(c(spntw$`higher level`[2], 
                                                              spntw$`higher level`[5],
                                                              spntw$`higher level`[3],
                                                              spntw$`higher level`[14],
                                                              spntw$`higher level`[20]), TRUE)
  outsp.site[n2+seq(nrow(spntw$`higher level`)),8]<-nestedcontr$`higher level` #other way around


  
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),1] <- as.character(sites[i])
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),2] <- try(as.character(rownames(spntw$`lower level`)), TRUE)
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),3:7] <- try(c(spntw$`lower level`[2], 
                                                                spntw$`lower level`[5],
                                                                spntw$`lower level`[3],
                                                                 spntw$`lower level`[14],
                                                                spntw$`lower level`[20]), TRUE)
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),8]<-nestedcontr$`lower level` #other way around

}  

str(out.site)
head(out.site)
out.site

#check correlation between both types of approaches in calculating corrected NODF
cor(out.site[-1,6:7], method = "spearman")

head(outsp.site)
str(outsp.site)

head(outsp.pl.site)
str(outsp.pl.site)

write.table(out.site[2:nrow(out.site),], file= "SITE_network_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.site[2:nrow(outsp.site),], file= "SITE_species_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.pl.site[2:nrow(outsp.pl.site),], file= "SITE_plant_species_level_metrics.csv", row.names= FALSE, sep= ",")


