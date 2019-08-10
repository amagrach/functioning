
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

##1st create a bipartite matrix for each site MERGING ALL ROUNDS TOGETHER and extract a set of metrics at the network level

sites <- unique(d$Site_ID)


out.site <- data.frame(Site_id = NA,  
                       species.poll=NA, species.pl=NA, 
                       functional.comp.poll=NA, functional.comp.plant=NA, 
                       nodf.song=NA)

outsp.site <- data.frame(Site_id = NA, species=NA,  norm_degree = NA,
                         weigh_closeness=NA, d=NA)


outsp.pl.site <- data.frame(Site_id = NA, species=NA,  norm_degree = NA,
                            weigh_closeness=NA, d=NA)



webs <- list()
for(i in 1:length(sites)){
  print(i)
  temp <- subset(d, Site_ID == sites[i])
  
  

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
  
   
  
  n <- nrow(out.site)
  webs[[n]] <- web4  
  n2 <- nrow(outsp.site)
  n3 <- nrow(outsp.pl.site)
  
  out.site[n + 1,1] <- as.character(sites[i])
  out.site[n + 1,2:5] <- c(ntw[20:21], ntw[42:43])
  out.site[n + 1,6] <- try(combined_NODF, TRUE)
 
  outsp.site[n2+seq(nrow(spntw$`higher level`)),1] <- as.character(sites[i])
  outsp.site[n2+seq(nrow(spntw$`higher level`)),2] <- try(as.character(rownames(spntw$`higher level`)), TRUE)
  outsp.site[n2+seq(nrow(spntw$`higher level`)),3:5] <- try(c(spntw$`higher level`[2], 
                                                              spntw$`higher level`[14],
                                                              spntw$`higher level`[20]), TRUE)


  
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),1] <- as.character(sites[i])
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),2] <- try(as.character(rownames(spntw$`lower level`)), TRUE)
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),3:5] <- try(c(spntw$`lower level`[2], 
                                                                 spntw$`lower level`[14],
                                                                spntw$`lower level`[20]), TRUE)

}  

str(out.site)
head(out.site)
out.site

head(outsp.site)
str(outsp.site)

head(outsp.pl.site)
str(outsp.pl.site)

write.table(out.site[2:nrow(out.site),], file= "SITE_network_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.site[2:nrow(outsp.site),], file= "SITE_species_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.pl.site[2:nrow(outsp.pl.site),], file= "SITE_plant_species_level_metrics.csv", row.names= FALSE, sep= ",")


