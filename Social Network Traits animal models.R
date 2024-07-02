
#### ANIMAL MODELS FOR SOCIAL NETWORK TRAITS #####

#REQUIRED LIBRARIES

library(dplyr)
library(asreml)
library(ggplot2)


# ---------- PREPARE DATA ----------------

#Read in metrics data and pedigree data

flockswknd <- readRDS("/ENTERYOURDIRECTORY/socialnetworkdatafull.rds")

flockswknd <- flockswknd %>% dplyr::mutate_at(c('id', 'nwkend', 'adjuv', 'year', 'sex', 'adjuv'), as.factor) 
flockswknd$id <- factor(toupper(levels(flockswknd$id))[flockswknd$id])

print(length(unique(flockswknd$id))) #1823

ped_data <- read.csv("/ENTERYOURDIRECTORY/GTITprunedpedigree.csv", na.strings=c("", "NA"))
ped_inv <- ainverse(ped_data)

# ---------- RUN MODELS ----------------

#TRAITS: mean flock size, degree, weighted degree, weighted eigencentrality, eigencentrality, betweenness


#2011

flockswknd2011 <- flockswknd %>% filter(year == 2011)
print(length(unique(flockswknd2011$id))) #1085 ids (9855 obs)

wknd2011fs <- asreml(fixed= mean.flock.size ~ 1, #Mean Group size
                     random=~id,
                     residual=~idv(units),
                     data=flockswknd2011,
                     na.action=na.method(x="omit", y="omit"))
summary(wknd2011fs)$varcomp
vpredict(wknd2011fs, repeatabilityfs ~ (V1)/(V1+V2)) #0.379 (0.013)

wknd2011deg <- asreml(fixed= degree ~ 1,      #Degree
                      random=~id,
                      residual=~idv(units),
                      data=flockswknd2011,
                      na.action=na.method(x="omit", y="omit"))
summary(wknd2011deg)$varcomp
vpredict(wknd2011deg, repeatabilitydeg ~ (V1)/(V1+V2)) #0.379 (0.0134)

wknd2011wdeg <- asreml(fixed= w.degree ~ 1,      #Weighted degree
                       random=~id,
                       residual=~idv(units),
                       data=flockswknd2011,
                       na.action=na.method(x="omit", y="omit"))
summary(wknd2011wdeg)$varcomp
vpredict(wknd2011wdeg, repeatabilitywdeg ~ (V1)/(V1+V2)) #0.375 (0.0134)

wknd2011wcentr <- asreml(fixed= w.eigen_cent ~ 1,      #Weighted centrality
                         random=~id,
                         residual=~idv(units),
                         data=flockswknd2011,
                         na.action=na.method(x="omit", y="omit"))
summary(wknd2011wcentr)$varcomp
vpredict(wknd2011wcentr, repeatabilitywcentr ~ (V1)/(V1+V2)) #0.086 (0.007)

wknd2011centr <- asreml(fixed= eigen_cent ~ 1,      #Centrality
                        random=~id,
                        residual=~idv(units),
                        data=flockswknd2011,
                        na.action=na.method(x="omit", y="omit"))
summary(wknd2011centr)$varcomp
vpredict(wknd2011centr, repeatabilitycentr ~ (V1)/(V1+V2)) #0.126 (0.009)

wknd2011betw <- asreml(fixed= betweenness ~ 1,     #Betweenness
                       random=~id,
                       residual=~idv(units),
                       data=flockswknd2011,
                       na.action=na.method(x="omit", y="omit"))
summary(wknd2011betw)$varcomp
vpredict(wknd2011betw, repeatabilitybetw ~ (V1)/(V1+V2)) #0.0551 (0.006)

#heritabilities

wknd2011fsped <- asreml(fixed= mean.flock.size ~ 1,           #Mean group size        
                        random=~ vm(id, ped_inv) + ide(id),    
                        residual=~idv(units),
                        data=flockswknd2011,
                        na.action=na.method(x="omit", y="omit"))
summary(wknd2011fsped)$varcomp

vpredict(wknd2011fsped, repeatabilityfs ~ (V1+V2)/(V1+V2+V3)) #0.379 (0.013)
vpredict(wknd2011fsped, h2fs ~ (V1)/(V1+V2+V3)) #4.421435e-08 (9.651839e-10)

wknd2011degped <- asreml(fixed= degree ~ 1,                   #Degree
                         random=~ vm(id, ped_inv) + ide(id),    
                         residual=~idv(units),
                         data=flockswknd2011,
                         na.action=na.method(x="omit", y="omit"))
summary(wknd2011degped)$varcomp

vpredict(wknd2011degped, repeatabilitydeg ~ (V1+V2)/(V1+V2+V3)) #0.379 (0.013)
vpredict(wknd2011degped, h2deg ~ (V1)/(V1+V2+V3)) #5.306944e-08 (1.164433e-09)

wknd2011wdegped <- asreml(fixed= w.degree ~ 1,                   #Weighted degree
                          random=~ vm(id, ped_inv) + ide(id),    
                          residual=~idv(units),
                          data=flockswknd2011,
                          na.action=na.method(x="omit", y="omit"))
summary(wknd2011wdegped)$varcomp

vpredict(wknd2011wdegped, repeatabilitywdeg ~ (V1+V2)/(V1+V2+V3)) #0.375 (0.013)
vpredict(wknd2011wdegped, h2wdeg ~ (V1)/(V1+V2+V3)) #5.114394e-08 (1.116158e-09)

wknd2011wcentrped <- asreml(fixed= w.eigen_cent ~ 1,                 #Weighted centrality  
                            random=~ vm(id, ped_inv) + ide(id),    
                            residual=~idv(units),
                            data=flockswknd2011,
                            na.action=na.method(x="omit", y="omit"))
summary(wknd2011wcentrped)$varcomp

vpredict(wknd2011wcentrped, repeatabilitywcentr ~ (V1+V2)/(V1+V2+V3)) #0.09 (0.008)
vpredict(wknd2011wcentrped, h2wcentr ~ (V1)/(V1+V2+V3)) #0.02 (0.01)

wknd2011centrped <- asreml(fixed= eigen_cent ~ 1,                  #Centrality 
                           random=~ vm(id, ped_inv) + ide(id),    
                           residual=~idv(units),
                           data=flockswknd2011,
                           na.action=na.method(x="omit", y="omit"))
summary(wknd2011centrped)$varcomp

vpredict(wknd2011centrped, repeatabilitycentr ~ (V1+V2)/(V1+V2+V3)) #0.131 (0.009)
vpredict(wknd2011centrped, h2centr ~ (V1)/(V1+V2+V3)) #0.028 (0.012)

wknd2011betwped <- asreml(fixed= betweenness ~ 1,                   #Betweenness
                          random=~ vm(id, ped_inv) + ide(id),    
                          residual=~idv(units),
                          data=flockswknd2011,
                          na.action=na.method(x="omit", y="omit"))
summary(wknd2011betwped)$varcomp

vpredict(wknd2011betwped, repeatabilitybetw ~ (V1+V2)/(V1+V2+V3)) #0.057 (0.007)
vpredict(wknd2011betwped, h2betw ~ (V1)/(V1+V2+V3)) #0.0123 (0.009)

#..............................................................................................

#2012

flockswknd2012 <- flockswknd %>% filter(year == 2012)
print(length(unique(flockswknd2012$id))) #720 birds (6538 obs)

#repeatabilities

wknd2012fs <- asreml(fixed= mean.flock.size ~ 1,
                     random=~id,
                     residual=~idv(units),
                     data=flockswknd2012,
                     na.action=na.method(x="omit", y="omit"))
summary(wknd2012fs)$varcomp
vpredict(wknd2012fs, repeatabilityfs ~ (V1)/(V1+V2)) #0.586 (0.015)

wknd2012deg <- asreml(fixed= degree ~ 1,
                      random=~id,
                      residual=~idv(units),
                      data=flockswknd2012,
                      na.action=na.method(x="omit", y="omit"))
summary(wknd2012deg)$varcomp
vpredict(wknd2012deg, repeatabilitydeg ~ (V1)/(V1+V2)) #0.574 (0.0155)

wknd2012wdeg <- asreml(fixed= w.degree ~ 1,
                       random=~id,
                       residual=~idv(units),
                       data=flockswknd2012,
                       na.action=na.method(x="omit", y="omit"))
summary(wknd2012wdeg)$varcomp
vpredict(wknd2012wdeg, repeatabilitywdeg ~ (V1)/(V1+V2)) #0.554 (0.0157)

wknd2012wcentr <- asreml(fixed= w.eigen_cent ~ 1,
                         random=~id,
                         residual=~idv(units),
                         data=flockswknd2012,
                         na.action=na.method(x="omit", y="omit"))
summary(wknd2012wcentr)$varcomp
vpredict(wknd2012wcentr, repeatabilitywcentr ~ (V1)/(V1+V2)) #0.259 (0.014)

wknd2012centr <- asreml(fixed= eigen_cent ~ 1,
                        random=~id,
                        residual=~idv(units),
                        data=flockswknd2012,
                        na.action=na.method(x="omit", y="omit"))
summary(wknd2012centr)$varcomp
vpredict(wknd2012centr, repeatabilitycentr ~ (V1)/(V1+V2)) #0.365 (0.016)

wknd2012betw <- asreml(fixed= betweenness ~ 1,
                       random=~id,
                       residual=~idv(units),
                       data=flockswknd2012,
                       na.action=na.method(x="omit", y="omit"))
summary(wknd2012betw)$varcomp
vpredict(wknd2012betw, repeatabilitybetw ~ (V1)/(V1+V2)) #0.0423 (0.008)

#heritabilities

wknd2012fsped <- asreml(fixed= mean.flock.size ~ 1,                   
                        random=~ vm(id, ped_inv) + ide(id),    
                        residual=~idv(units),
                        data=flockswknd2012,
                        na.action=na.method(x="omit", y="omit"))
summary(wknd2012fsped)$varcomp

vpredict(wknd2012fsped, repeatabilityfs ~ (V1+V2)/(V1+V2+V3)) #0.591 (0.015)
vpredict(wknd2012fsped, h2fs ~ (V1)/(V1+V2+V3)) #0.062 (0.051) 

wknd2012degped <- asreml(fixed= degree ~ 1,                   
                         random=~ vm(id, ped_inv) + ide(id),    
                         residual=~idv(units),
                         data=flockswknd2012,
                         na.action=na.method(x="omit", y="omit"))
summary(wknd2012degped)$varcomp

vpredict(wknd2012degped, repeatabilitydeg ~ (V1+V2)/(V1+V2+V3)) #0.577 (0.016)
vpredict(wknd2012degped, h2deg ~ (V1)/(V1+V2+V3)) #0.0499 (0.0491)

wknd2012wdegped <- asreml(fixed= w.degree ~ 1,                   
                          random=~ vm(id, ped_inv) + ide(id),    
                          residual=~idv(units),
                          data=flockswknd2012,
                          na.action=na.method(x="omit", y="omit"))
summary(wknd2012wdegped)$varcomp

vpredict(wknd2012wdegped, repeatabilitywdeg ~ (V1+V2)/(V1+V2+V3)) #0.557 (0.0136)
vpredict(wknd2012wdegped, h2wdeg ~ (V1)/(V1+V2+V3)) #0.043 (0.049)

wknd2012wcentrped <- asreml(fixed= w.eigen_cent ~ 1,                   
                            random=~ vm(id, ped_inv) + ide(id),    
                            residual=~idv(units),
                            data=flockswknd2012,
                            na.action=na.method(x="omit", y="omit"))
summary(wknd2012wcentrped)$varcomp

vpredict(wknd2012wcentrped, repeatabilitywcentr ~ (V1+V2)/(V1+V2+V3)) #0.265 (0.015)
vpredict(wknd2012wcentrped, h2wcentr ~ (V1)/(V1+V2+V3)) #0.048 (0.029)

wknd2012centrped <- asreml(fixed= eigen_cent ~ 1,                   
                           random=~ vm(id, ped_inv) + ide(id),    
                           residual=~idv(units),
                           data=flockswknd2012,
                           na.action=na.method(x="omit", y="omit"))
summary(wknd2012centrped)$varcomp

vpredict(wknd2012centrped, repeatabilitycentr ~ (V1+V2)/(V1+V2+V3)) #0.374 (0.017)
vpredict(wknd2012centrped, h2centr ~ (V1)/(V1+V2+V3)) #0.103 (0.035)

wknd2012centrped$loglik #7063.051
wknd2012centr$loglik # 7057.361
1- pchisq(2 * (wknd2012centrped$loglik - wknd2012centr$loglik),1) #0.0007 p value

wknd2012betwped <- asreml(fixed= betweenness ~ 1,                   
                          random=~ vm(id, ped_inv) + ide(id),    
                          residual=~idv(units),
                          data=flockswknd2012,
                          na.action=na.method(x="omit", y="omit"))
summary(wknd2012betwped)$varcomp

vpredict(wknd2012betwped, repeatabilitybetw ~ (V1+V2)/(V1+V2+V3)) #0.042 (0.008)
vpredict(wknd2012betwped, h2betw ~ (V1)/(V1+V2+V3)) #1.683863e-08 (2.976612e-10)

#..............................................................................................

#2013

flockswknd2013 <- flockswknd %>% filter(year == 2013)
print(length(unique(flockswknd2013$id))) #797 birds (7348 obs)

#repeatabilities

wknd2013fs <- asreml(fixed= mean.flock.size ~ 1,
                     random=~id,
                     residual=~idv(units),
                     data=flockswknd2013,
                     na.action=na.method(x="omit", y="omit"))
summary(wknd2013fs)$varcomp
vpredict(wknd2013fs, repeatabilityfs ~ (V1)/(V1+V2)) #0.545 (0.015)

wknd2013deg <- asreml(fixed= degree ~ 1,
                      random=~id,
                      residual=~idv(units),
                      data=flockswknd2013,
                      na.action=na.method(x="omit", y="omit"))
summary(wknd2013deg)$varcomp
vpredict(wknd2013deg, repeatabilitydeg ~ (V1)/(V1+V2)) #0.496 (0.0155)

wknd2013wdeg <- asreml(fixed= w.degree ~ 1,
                       random=~id,
                       residual=~idv(units),
                       data=flockswknd2013,
                       na.action=na.method(x="omit", y="omit"))
summary(wknd2013wdeg)$varcomp
vpredict(wknd2013wdeg, repeatabilitywdeg ~ (V1)/(V1+V2)) #0.537 (0.0151)

wknd2013wcentr <- asreml(fixed= w.eigen_cent ~ 1,
                         random=~id,
                         residual=~idv(units),
                         data=flockswknd2013,
                         na.action=na.method(x="omit", y="omit"))
summary(wknd2013wcentr)$varcomp
vpredict(wknd2013wcentr, repeatabilitywcentr ~ (V1)/(V1+V2)) #0.046 (0.007)

wknd2013centr <- asreml(fixed= eigen_cent ~ 1,
                        random=~id,
                        residual=~idv(units),
                        data=flockswknd2013,
                        na.action=na.method(x="omit", y="omit"))
summary(wknd2013centr)$varcomp
vpredict(wknd2013centr, repeatabilitycentr ~ (V1)/(V1+V2)) #0.135 (0.01)

wknd2013betw <- asreml(fixed= betweenness ~ 1,
                       random=~id,
                       residual=~idv(units),
                       data=flockswknd2013,
                       na.action=na.method(x="omit", y="omit"))
summary(wknd2013betw)$varcomp
vpredict(wknd2013betw, repeatabilitybetw ~ (V1)/(V1+V2)) #0.076 (0.008)

#heritabilities

wknd2013fsped <- asreml(fixed= mean.flock.size ~ 1,                   
                        random=~ vm(id, ped_inv) + ide(id),    
                        residual=~idv(units),
                        data=flockswknd2013,
                        na.action=na.method(x="omit", y="omit"))
summary(wknd2013fsped)$varcomp

vpredict(wknd2013fsped, repeatabilityfs ~ (V1+V2)/(V1+V2+V3)) #0.549 (0.015)
vpredict(wknd2013fsped, h2fs ~ (V1)/(V1+V2+V3)) #0.062 (0.041) 

wknd2013degped <- asreml(fixed= degree ~ 1,                   
                         random=~ vm(id, ped_inv) + ide(id),    
                         residual=~idv(units),
                         data=flockswknd2013,
                         na.action=na.method(x="omit", y="omit"))
summary(wknd2013degped)$varcomp

vpredict(wknd2013degped, repeatabilitydeg ~ (V1+V2)/(V1+V2+V3)) #0.499 (0.016)
vpredict(wknd2013degped, h2deg ~ (V1)/(V1+V2+V3)) #0.042 (0.036)

wknd2013wdegped <- asreml(fixed= w.degree ~ 1,                   
                          random=~ vm(id, ped_inv) + ide(id),    
                          residual=~idv(units),
                          data=flockswknd2013,
                          na.action=na.method(x="omit", y="omit"))
summary(wknd2013wdegped)$varcomp

vpredict(wknd2013wdegped, repeatabilitywdeg ~ (V1+V2)/(V1+V2+V3)) #0.546 (0.016)
vpredict(wknd2013wdegped, h2wdeg ~ (V1)/(V1+V2+V3)) #0.114 (0.042)

wknd2013wcentrped <- asreml(fixed= w.eigen_cent ~ 1,                   
                            random=~ vm(id, ped_inv) + ide(id),    
                            residual=~idv(units),
                            data=flockswknd2013,
                            na.action=na.method(x="omit", y="omit"))
summary(wknd2013wcentrped)$varcomp

vpredict(wknd2013wcentrped, repeatabilitywcentr ~ (V1+V2)/(V1+V2+V3)) #0.054 (0.008)
vpredict(wknd2013wcentrped, h2wcentr ~ (V1)/(V1+V2+V3)) #0.043 (0.012)

wknd2013centrped <- asreml(fixed= eigen_cent ~ 1,                   
                           random=~ vm(id, ped_inv) + ide(id),    
                           residual=~idv(units),
                           data=flockswknd2013,
                           na.action=na.method(x="omit", y="omit"))
summary(wknd2013centrped)$varcomp

vpredict(wknd2013centrped, repeatabilitycentr ~ (V1+V2)/(V1+V2+V3)) #0.141 (0.011)
vpredict(wknd2013centrped, h2centr ~ (V1)/(V1+V2+V3)) #0.037 (0.017)

wknd2013betwped <- asreml(fixed= betweenness ~ 1,                   
                          random=~ vm(id, ped_inv) + ide(id),    
                          residual=~idv(units),
                          data=flockswknd2013,
                          na.action=na.method(x="omit", y="omit"))
summary(wknd2013betwped)$varcomp

vpredict(wknd2013betwped, repeatabilitybetw ~ (V1+V2)/(V1+V2+V3)) #0.077 (0.009)
vpredict(wknd2013betwped, h2betw ~ (V1)/(V1+V2+V3)) #0.007 (0.012)

#..............................................................................................

#between year (year as fixed effect)

#repeatabilities

allfs <- asreml(fixed= mean.flock.size ~ 1 + year,
                random=~id,
                residual=~idv(units),
                data=flockswknd,
                na.action=na.method(x="omit", y="omit"))
summary(allfs)$varcomp
vpredict(allfs, repeatabilityfs ~ (V1)/(V1+V2)) #0.4236565 (0.01006839)

alldeg <- asreml(fixed= degree ~ 1 + year,
                 random=~id,
                 residual=~idv(units),
                 data=flockswknd,
                 na.action=na.method(x="omit", y="omit"))
summary(alldeg)$varcomp
vpredict(alldeg, repeatabilitydeg~ (V1)/(V1+V2)) #0.4140583 (0.01012216)

allwdeg <- asreml(fixed= w.degree ~ 1 + year,
                  random=~id,
                  residual=~idv(units),
                  data=flockswknd,
                  na.action=na.method(x="omit", y="omit"))
summary(allwdeg)$varcomp
vpredict(allwdeg, repeatabilitywdeg ~ (V1)/(V1+V2)) #0.416 (0.01)

allwcentr <- asreml(fixed= w.eigen_cent ~ 1 + year,
                    random=~id,
                    residual=~idv(units),
                    data=flockswknd,
                    na.action=na.method(x="omit", y="omit"))
summary(allwcentr)$varcomp
vpredict(allwcentr, repeatabilitywcentr ~ (V1)/(V1+V2)) #0.097 (0.005)

allcentr <- asreml(fixed= eigen_cent ~ 1 + year,
                   random=~id,
                   residual=~idv(units),
                   data=flockswknd,
                   na.action=na.method(x="omit", y="omit"))
summary(allcentr)$varcomp
vpredict(allcentr, repeatabilitycentr ~ (V1)/(V1+V2)) #0.169 (0.007)

allbetw <- asreml(fixed= betweenness ~ 1 + year,
                  random=~id,
                  residual=~idv(units),
                  data=flockswknd,
                  na.action=na.method(x="omit", y="omit"))
summary(allbetw)$varcomp
vpredict(allbetw, repeatabilitybetw ~ (V1)/(V1+V2)) #0.049 (0.004)

#heritabilities

allfsped <- asreml(fixed= mean.flock.size ~ 1 + year,                   
                   random=~ vm(id, ped_inv) + ide(id),    
                   residual=~idv(units),
                   data=flockswknd,
                   na.action=na.method(x="omit", y="omit"))
summary(allfsped)$varcomp
vpredict(allfsped, repeatabilityfs ~ (V1+V2)/(V1+V2+V3)) #0.4249323 (0.01032617)
vpredict(allfsped, h2fs ~ (V1)/(V1+V2+V3)) #0.01082405 (0.01704833)

alldegped <- asreml(fixed= degree ~ 1 + year,                   
                    random=~ vm(id, ped_inv) + ide(id),    
                    residual=~idv(units),
                    data=flockswknd,
                    na.action=na.method(x="omit", y="omit"))
summary(alldegped)$varcomp
vpredict(alldegped, repeatabilityfs ~ (V1+V2)/(V1+V2+V3)) #0.4140584 (0.01012216)
vpredict(alldegped, h2fs ~ (V1)/(V1+V2+V3)) #6.075287e-08 (1.049013e-09)

allwdegped <- asreml(fixed= w.degree ~ 1 + year,                   
                     random=~ vm(id, ped_inv) + ide(id),    
                     residual=~idv(units),
                     data=flockswknd,
                     na.action=na.method(x="omit", y="omit"))
summary(allwdegped)$varcomp
vpredict(allwdegped, repeatabilitywdeg ~ (V1+V2)/(V1+V2+V3)) #0.418 (0.010)
vpredict(allwdegped, h2wdeg ~ (V1)/(V1+V2+V3)) #0.012 (0.017)

allwcentr <- asreml(fixed= w.eigen_cent ~ 1 + year,                   
                    random=~ vm(id, ped_inv) + ide(id),    
                    residual=~idv(units),
                    data=flockswknd,
                    na.action=na.method(x="omit", y="omit"))
summary(allwcentr)$varcomp
vpredict(allwcentr, repeatabilitywcentr ~ (V1+V2)/(V1+V2+V3)) #0.099 (0.005)
vpredict(allwcentr, h2wcentr ~ (V1)/(V1+V2+V3)) #0.010 (0.006)

allcentr <- asreml(fixed= eigen_cent ~ 1 + year,                   
                   random=~ vm(id, ped_inv) + ide(id),    
                   residual=~idv(units),
                   data=flockswknd,
                   na.action=na.method(x="omit", y="omit"))
summary(allcentr)$varcomp
vpredict(allcentr, repeatabilitycentr ~ (V1+V2)/(V1+V2+V3)) #0.172 (0.007)
vpredict(allcentr, h2centr ~ (V1)/(V1+V2+V3)) #0.022 (0.009)

allbetw <- asreml(fixed= betweenness ~ 1 + year,                   
                  random=~ vm(id, ped_inv) + ide(id),    
                  residual=~idv(units),
                  data=flockswknd,
                  na.action=na.method(x="omit", y="omit"))
summary(allbetw)$varcomp
vpredict(allbetw, repeatabilitybetw ~ (V1+V2)/(V1+V2+V3)) #0.050 (0.004)
vpredict(allbetw, h2betw ~ (V1)/(V1+V2+V3)) #0.008 (0.004)

#------------------------------------------------------------------------------------

# ---------- FIGURES ----------------

#FIGURE 4: Repeatability and Heritability estimates from Social Network Trait models

#Read in file with saved estimates from models 
heritsocnet <- read.csv("/ENTERYOURDIRECTORY/socnet model estimates.csv", na.strings=c("", "NA"))
herit2011 <- heritsocnet %>% filter(Year == "2011")

#Following is the plot for 2011 models. Replace the 2011 estimates with other years to obtain plots for the corresponding winters

socnet2011 <- ggplot(herit2011, aes(x = Estimate, y = Trait, color = EstimateType)) +
  geom_point(size=8) +
  geom_errorbarh(aes(xmin = Estimate - SE, xmax = Estimate + SE),
                 position = position_dodge(width = 1.0), width = 0.5, color="black", height=0.1) +
  labs(x = "Estimate", y = "Trait", color = "Estimate Type") +
  theme(legend.position = "none") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))+
  scale_color_manual(values = c("#d96459", "#1FBBC6")) + xlim(0, 0.6) +
  theme(axis.text = element_text(size = 25),
        axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title.x = element_text(size = 25, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        plot.title = element_text(size = 25, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        panel.grid.major.y = element_line(color = "grey", linetype = "dotted"))  # Add horizontal grid lines
+
  ggtitle("2011 winter")
