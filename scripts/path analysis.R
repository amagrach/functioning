

library(lavaan)
library(semPlot)
library(OpenMx)
library(tidyverse)
library(knitr)
library(kableExtra)
library(GGally)

# Organizing package information for table
packages <- c("tidyverse", "knitr", "kableExtra", "lavaan", "semPlot", "OpenMx", "GGally")
display <- c("Package","Title", "Maintainer", "Version", "URL")
table <- matrix(NA, 1, NROW(display), dimnames = list(1, display))
for(i in 1:NROW(packages)){
  list <- packageDescription(packages[i])
  table <- rbind(table, matrix(unlist(list[c(display)]), 1, NROW(display), byrow = T))
}
table[,NROW(display)] <- stringr::str_extract(table[,NROW(display)], ".+,")

# Table of packages
kable(table[-1,], format = "html", align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))




#species level 


model <-'
tot.fruitset ~ norm_degree.sc  + tot.visits.sc + 
weigh_closeness.sc + niche.overlap.sc
'


fit <- cfa(model, data = f2_scaled)
summary(fit, fit.measures = TRUE, standardized=T,rsquare=T)

semPaths(fit, 'std', layout = 'circle')

semPaths(fit,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)



ggcorr(f2_scaled[,c(14, 26:29)], nbreaks = 6, label = T, low = "red3", high = "green3", 
       label_round = 2, name = "Correlation Scale", label_alpha = T, hjust = 0.75) +
  ggtitle(label = "Correlation Plot") +
  theme(plot.title = element_text(hjust = 0.6))





#community level

model <-'
mean.fruitset ~ species.poll.sc + tot.visits.sc + nodf.song.sc + 
functional.comp.poll.sc
'


fit <- cfa(model, data = f.agg)
summary(fit, fit.measures = TRUE, standardized=T,rsquare=T)

semPaths(fit, 'std', layout = 'circle')

semPaths(fit,"std",layout = 'tree', edge.label.cex=1.5, curvePivot = TRUE)





#equity


f6red <- f6[, c(1:2, 7:10)]
colnames(f6red) <- c("Site", "Eq", "N", "F", "D", "V")

model <-'
Eq ~ D + V + N + F'


fit <- cfa(model, data = f6red)
summary(fit, fit.measures = TRUE, standardized=T,rsquare=T)

semPaths(fit, 'std', layout = 'circle')

semPaths(fit,"std",layout = 'tree', edge.label.cex=1.5, curvePivot = TRUE)
