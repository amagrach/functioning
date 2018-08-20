
#First attempt to load the data, etc...

#prepare stuff-----
#If you don't have devtools install it
#install.packages("devtools")
library(devtools)
#Then install the data package. Reprot issues 
install_github("ibartomeus/BeeFunData", ref = "fun")
library(BeeFunData)

#make the data available
data(all_interactions)
head(all_interactions)
#See metadata
?all_interactions
data(sites)
sites #lat longs
#load plant functioning data
data(fruitset)
?fruitset
head(fruitset)
data(seedset)
?seedset
data(single_visits)
?single_visits

#analysis...----

setwd("~/Dropbox/Beefun/Beefun/functioning")

d<-all_interactions
f<-fruitset
s<-seedset
si<-single_visits


##1st create a bipartite matrix for each site MERGING ALL ROUNDS TOGETHER and extract a set of metrics at the network level

#############################################################################
#####################Network and SPECIES LEVEL NETWORK METRICS GLOBAL NETWORK PER SITE##############################
#############################################################################

## super importante cargar primero igraph antes que network porque se sobrescriben partes de las librerias
## y el script no funcionaria bien.

require(igraph)
require(network)
require(sna)
require(ggplot2)
require(sna)
require(ergm)
require(bipartite)
require(dplyr)
require(reshape)
require(wesanderson)
# require(ESM)


sites <- unique(d$Site_ID)


out.site <- data.frame(Site_id = NA, links_sps = NA, num_compartments = NA,weighted_NODF  = NA,
                       link_dens=NA, weigh_connectance=NA, Shannon=NA,
                       int_evenness=NA,  H2=NA, 
                       species.poll=NA, species.pl=NA, 
                       niche.overlap.poll=NA, niche.overlap.plant=NA, robustness.poll=NA, robustness.plant=NA,
                       functional.comp.poll=NA, functional.comp.plant=NA,
                       modularity=NA, z=NA, links.corr=NA, nodf.corr=NA, linkdens.corr=NA, conn.corr=NA,
                       even.corr=NA, h2.corr=NA)

outsp.site <- data.frame(Site_id = NA, species=NA,  norm_degree = NA,strength= NA, nestedra= NA,weigh_betweeness  = NA, 
                         weigh_closeness=NA, d=NA, c=NA, z=NA)#, specificity=NA)


outsp.pl.site <- data.frame(Site_id = NA, species=NA,  norm_degree = NA,strength= NA, nestedra= NA,
                            weigh_betweeness  = NA, 
                            weigh_closeness=NA, d=NA, c=NA, z=NA)#, specificity=NA)


webs <- list()
for(i in 1:length(sites)){
  temp <- subset(d, Site_ID == sites[i])
  # for(j in 1:length(rounds)){
  #   temp2 <- subset(temp, Round == rounds[j])
  #   if (nrow(temp2) == 0) next
  
  web <- table(temp$Plant_gen_sp, temp$Pollinator_gen_sp)
  web2<-as.data.frame.array(web)
  #p1<-plotweb(web)
  spntw <- try(specieslevel(web), TRUE)
  ntw <- networklevel(web)
  #if there is an error skip this round and continue
  if(isTRUE(class(spntw)=="try-error")) {next} 
  
  #para calcular modularity de la red
  mod<-try(computeModules(web2), TRUE)
  
  #null model 
  nulls <- try(nullmodel(web, N=100, method="r2d"), TRUE)
  ntw.null <- try(sapply(nulls, networklevel), TRUE)
  
  #null value of modularity
  modules.nulls <- try(sapply(nulls, computeModules), TRUE)
  like.nulls <- try(sapply(modules.nulls, function(x) x@likelihood), TRUE)
  z <- try((mod@likelihood - mean(like.nulls))/sd(like.nulls), TRUE)
  
  #species-level c and z values for plants and poll
  
  cz<-try(czvalues(mod, level="higher", weighted=TRUE), TRUE)
  cz.pl<-try(czvalues(mod, level="lower", weighted=TRUE), TRUE)
  
  # spec<-getspe(web)
  
  
  n <- nrow(out.site)
  webs[[n]] <- web
  n2 <- nrow(outsp.site)
  n3 <- nrow(outsp.pl.site)
  
  out.site[n + 1,1] <- as.character(sites[i])
  # out[n + 1,2] <- as.character(rounds[j])
  out.site[n + 1,2:17] <- c(ntw[3:4],ntw[10], ntw[13:14], ntw[16:17],
                            ntw[19:21], ntw[28:29], ntw[40:43])
  out.site[n + 1,18] <- try(mod@likelihood, TRUE)
  out.site[n + 1,19] <- try(z, TRUE)
  out.site[n + 1,20] <- try((ntw[3] - mean(ntw.null[3,1:100])/sd(ntw.null[3,1:100])), TRUE)
  out.site[n + 1,21] <- try((ntw[9] - mean(ntw.null[9,1:100])/sd(ntw.null[9,1:100])), TRUE)
  out.site[n + 1,22] <- try((ntw[12] - mean(ntw.null[12,1:100])/sd(ntw.null[12,1:100])), TRUE)
  out.site[n + 1,23] <- try((ntw[13] - mean(ntw.null[13,1:100])/sd(ntw.null[13,1:100])), TRUE)
  out.site[n + 1,24] <- try((ntw[16] - mean(ntw.null[15,1:100])/sd(ntw.null[15,1:100])), TRUE)
  out.site[n + 1,25] <- try((ntw[18] - mean(ntw.null[16,1:100])/sd(ntw.null[16,1:100])), TRUE)
  
  
  
  outsp.site[n2+seq(nrow(spntw$`higher level`)),1] <- as.character(sites[i])
  # outsp[n2+seq(nrow(spntw$`higher level`)),2] <- as.character(rounds[j])
  outsp.site[n2+seq(nrow(spntw$`higher level`)),2] <- try(as.character(rownames(spntw$`higher level`)), TRUE)
  outsp.site[n2+seq(nrow(spntw$`higher level`)),3:8] <- try(c(spntw$`higher level`[2:3], 
                                                              spntw$`higher level`[5],
                                                              spntw$`higher level`[12], spntw$`higher level`[14],
                                                              spntw$`higher level`[20]), TRUE)
  outsp.site[n2+seq(nrow(spntw$`higher level`)),9:10] <- try(c(cz$c,
                                                               cz$z), TRUE)
  # outsp[n2+seq(nrow(web)),12] <- spec #calculates specificity
  
  
  
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),1] <- as.character(sites[i])
  # outsp.pl[n2+seq(nrow(spntw$`lower level`)),2] <- as.character(rounds[j])
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),2] <- try(as.character(rownames(spntw$`lower level`)), TRUE)
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),3:8] <- try(c(spntw$`lower level`[2:3], 
                                                                spntw$`lower level`[5],
                                                                spntw$`lower level`[12], spntw$`lower level`[14],
                                                                spntw$`lower level`[20]), TRUE)
  outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),9:10] <- try(c(cz.pl$c,
                                                                 cz.pl$z), TRUE)
  # outsp.pl[n2+seq(nrow(web)),12] <- spec #calculates specificity
  
}  
# }  

str(out.site)
head(out.site)
out.site

head(outsp.site)
str(outsp.site)

head(outsp.pl.site)
str(outsp.pl.site)

write.table(out.site[2:17,], file= "SITE_network_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.site[2:763,], file= "SITE_species_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.pl.site[2:209,], file= "SITE_plant_species_level_metrics.csv", row.names= FALSE, sep= ",")

