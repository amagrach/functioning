#calculate average niche overlap

library(reshape2)
library(spaa)
sites <- unique(d$Site_ID)

out.niche <- data.frame(Site_id = NA, row = NA, col = NA,value  = NA)

webs <- list()
for(i in 1:length(sites)){
  temp <- subset(d, Site_ID == sites[i])
  temp2<-droplevels(temp)
 
  web <- table(temp2$Pollinator_gen_sp, temp2$Plant_gen_sp)

  ni<-niche.overlap(web, method="morisita")
  df <- melt(as.matrix(ni), varnames = c("row", "col"))
  df$row<-as.character(df$row)
  df$col<-as.character(df$col)
  
  n <- nrow(out.niche)
 
  out.niche[n + seq(nrow(df)),1] <- as.character(sites[i])
  out.niche[n + seq(nrow(df)),2:3] <- (df[,1:2])
  out.niche[n + seq(nrow(df)),4] <-df[,3]
  
}



out.niche.agg <- ddply(out.niche, c("Site_id", "col"), summarise,
               mean.niche.over    = mean(value))
