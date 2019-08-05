library(BeeFunData)
d<-all_interactions
library(plyr)
library(dplyr)


  
  sites <- unique(d$Site_ID)


out <- data.frame(Site_id = NA, plant.sps = NA, poll.sps = NA)

  
   for(i in 1:length(sites)){
    temp <- subset(d, Site_ID == sites[i])
    temp2<-droplevels(temp)
    pl<-nlevels(unique(temp2$Plant_gen_sp))
    poll<-nlevels(unique(temp2$Pollinator_gen_sp))
    
    n <- nrow(out)
    
    out[n + 1,1] <- as.character(sites[i])
    out[n + 1,2] <- pl
    out[n + 1,3] <- poll
   
   }

out


out2<-out[-1,]

cor(out2$plant.sps, out2$poll.sps, method = ("pearson"))

