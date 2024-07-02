
#### ANIMAL MODELS WITH NATAL EFFECTS #####

#REQUIRED LIBRARIES

library(dplyr)
library(asreml)
library(ggplot2)
library(tidyr)
library(sf)
library(Matrix)


# ---------- PREPARE DATA ----------------

#MODELS WERE CREATED FOR GROUPSIZE AND DEGREE 
#Focal datasets:

natalflocks <- readRDS("/ENTERYOURDIRECTORY/flockingevents_natalinfo.rds")
natalwkndflocks <- readRDS("/ENTERYOURDIRECTORY/socialnetwork_natalinfo.rds")

#Creating matrices for spatial proximity and natal environmental similarity

natalhabitatdata <- read.csv(
  file.path("/ENTERYOURDIRECTORY/natalhabitatdata.csv"),na.strings = c("", "NA"))


natalhabitatdata <- natalhabitatdata %>% 
  mutate_at(
    c("id", "nestbox"),
    as.factor
  )

#NATAL ENVIRONMENTAL SIMILARITY MATRIX

#getting relevant columns
env_var <- natalhabitatdata[, c(
  "altitude_m", "northness", "edge_edi",
  "no_trees_75m", "area_polygon_sqrt"
)]
# Centering the data, and scaling by standard deviation
env_var <- scale(env_var, center = TRUE, scale = T)
rownames(env_var) <- natalhabitatdata$id

#calculating the euclidean distance between each individual's environment parameters
env_var_euc <- as.matrix(dist(env_var,
                              method = "euclidean",
                              diag = TRUE,
                              upper = T
))

#Scaling the values between 0 and 1
env_euc_scal <- 1 - env_var_euc / max(env_var_euc, na.rm = TRUE)
colnames(env_euc_scal) <-
  rownames(env_euc_scal) <-
  natalhabitatdata$id

#env_euc_scal is the raw matrix
#making positive definite so it can be run in animal model
env_euc_scal_PD <- Matrix::nearPD(env_euc_scal)
natalenv_matrix <- env_euc_scal_PD$mat


#SPATIAL MATRIX 

# making coordinates into geometries
spatial_geom <- sf::st_as_sf(natalhabitatdata,
                             coords = c("x", "y"),
                             crs = 27700
)
# keeping just geometry
geom <- spatial_geom[, c("geometry")]
rownames(geom) <- spatial_geom$id
# making distance matrix
spat_mat <- as.matrix(sf::st_distance(geom, geom))
colnames(spat_mat) <- rownames(spat_mat) <- spatial_geom$id

#make it to a matrix and array
spat_mat_conv <- as.numeric(spat_mat)
dim(spat_mat_conv) <- c(
  nrow(spatial_geom),
  nrow(spatial_geom)
) 
colnames(spat_mat_conv) <- rownames(spat_mat_conv) <- spatial_geom$id

#spat_mat_conv is the raw spatial matrix

#scaling the matrix - diagonal is 1
spat_mat_conv_sc <- 1 - spat_mat_conv / max(spat_mat_conv, na.rm = TRUE)
spat_mat_conv_sc_PD <- Matrix::nearPD(spat_mat_conv_sc)
natalspat_matrix <- spat_mat_conv_sc_PD$mat

#Pedigree data for animal models:

ped_data <- read.csv("/ENTERYOURDIRECTORY/GTITprunedpedigree.csv", na.strings=c("", "NA"))
ped_inv <- ainverse(ped_data) #inverse of pedigree data - needed for ASReml


# ---------- RUN MODELS ----------------

#### ------- Groupsize natal models - all 3 winters with year as fixed effect ------

natalflocks <- natalflocks %>% dplyr::mutate_at(c('id', 'sex', 'adjuv', 'born', 'logger', 'year', 'nwkend', 'flock'), as.factor) 
print(length(unique(natalflocks$id))) #938

#only id (Model a)
natalbaseid <- asreml(fixed= flock.size ~ year,  #### only id as random effect
                      random=~id,
                      residual=~idv(units),
                      data=natalflocks,
                      na.action=na.method(x="omit", y="omit"), workspace = 64e6)
summary(natalbaseid)$varcomp
vpredict(natalbaseid, repeatabilityid ~ V1/(V1+V2)) #0.206 (0.008)

#heritability model (Model b)
natalallheritonly <- asreml(fixed= flock.size ~ year,
                            random=~ vm(id, ped_inv) + ide(id),    ####
                            residual=~idv(units),
                            data=natalflocks,
                            na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                            )
summary(natalallheritonly)$varcomp
vpredict(natalallheritonly, herit ~ V1/(V1+V2+V3)) #0.031 (0.014)
vpredict(natalallheritonly, id ~ (V1+V2)/(V1+V2+V3)) #0.207 (0.008)

#natal section (Model c)
natalsectall <- asreml(fixed= flock.size ~ year,  #### id and section as re
                       random=~id + Section,
                       residual=~idv(units),
                       data=natalflocks,
                       na.action=na.method(x="omit", y="omit"), workspace = 64e6)

summary(natalsectall)$varcomp
vpredict(natalsectall, sect ~ V1/(V1+V2+V3)) #0.008 (0.005)
vpredict(natalsectall, id ~ V2/(V1+V2+V3)) #0.200 (0.008)

#brood identity (Model d)
natalpnumall <- asreml(fixed= flock.size ~ year,  #### id and pnum as re
                       random=~id + Pnum,
                       residual=~idv(units),
                       data=natalflocks,
                       na.action=na.method(x="omit", y="omit"), workspace = 64e6)
summary(natalpnumall)$varcomp
vpredict(natalpnumall, pnum ~ V1/(V1+V2+V3)) #0.022 (0.010)
vpredict(natalpnumall, id ~ V2/(V1+V2+V3)) #0.185 (0.008)

#both section and brood identity (Model e)
natalbothall <- asreml(fixed= flock.size ~ year,  #### id section, pnum as re
                       random=~id + Pnum + Section,
                       residual=~idv(units),
                       data=natalflocks,
                       na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
                       )
summary(natalbothall)$varcomp
vpredict(natalbothall, sect ~ V1/(V1+V2+V3+V4)) #0.008 (0.005)
vpredict(natalbothall, pnum ~ V2/(V1+V2+V3+V4)) #0.015 (0.010)
vpredict(natalbothall, id ~ V3/(V1+V2+V3+V4)) #0.186 (0.012)


#herit + both natal effects + logger (Model f)
natalheritbothall <- asreml(fixed= flock.size ~ year,
                            random=~ vm(id, ped_inv) + ide(id) + Section + Pnum + logger,    ####
                            residual=~idv(units),
                            data=natalflocks,
                            na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                            )
summary(natalheritbothall)$varcomp
vpredict(natalheritbothall, sect ~ (V1)/(V1+V2+V3+V4+V5+V6)) #0.0001 (0.0006)
vpredict(natalheritbothall, logger ~ (V2)/(V1+V2+V3+V4+V5+V6)) #0.260 (0.034)
vpredict(natalheritbothall, pnum ~ (V3)/(V1+V2+V3+V4+V5+V6)) #0.006 (0.005)
vpredict(natalheritbothall, herit ~ (V4)/(V1+V2+V3+V4+V5+V6)) #0.0006 (0.007)
vpredict(natalheritbothall, id ~ (V4+V5)/(V1+V2+V3+V4+V5+V6)) #0.076 (0.007)


##Testing significance of random effects
natalsectall$loglik #-8076.129
natalpnumall$loglik #-8077.173
natalbothall$loglik
natalbaseid$loglik

#sect - id
1- pchisq(2 * (natalsectall$loglik - natalbaseid$loglik),1) #0 natal section improves model significantly
#pnum - id
1- pchisq(2 * (natalpnumall$loglik - natalbaseid$loglik),1) #0 brood identity improves model significantly

#both - sect for pnum
1- pchisq(2 * (natalbothall$loglik - natalsectall$loglik),1) #0.148 Pnum not significant
1- pchisq(2 * (natalbothall$loglik - natalpnumall$loglik),1) #0 sect significant

1- pchisq(2 * (natalbothall$loglik - natalbaseid$loglik),2) #0 both effects significantly improve the model

#for models with matrices

# create another column to associate with second matrix

natalflocks$id2 <- natalflocks$id
natalflocks <- natalflocks %>% dplyr::mutate_at(c('id', 'sex', 'adjuv', 'born', 'logger', 'year', 'nwkend', 'flock', 'id2'), as.factor) 

#pedigree, id and spatial matrix (Model g)
fsheritspatialmatrix <- asreml(fixed= flock.size ~ year,
                               random=~ vm(id, ped_inv) + ide(id) + vm(id2, natalspat_matrix),    ####
                               residual=~idv(units),
                               data=natalflocks,
                               na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                              )

summary(fsheritspatialmatrix)$varcomp
vpredict(fsheritspatialmatrix, spatialdistance ~ V1/(V1+V2+V3+V4)) #0.02717383 0.01611464
vpredict(fsheritspatialmatrix, id ~ V3/(V1+V2+V3+V4)) # 0.191272 0.01294572
vpredict(fsheritspatialmatrix, herit ~ V2/(V1+V2+V3+V4)) #0.0004429563 0.01083575

#pedigree, id and natal env similarity matrix (Model h)
fsheritenvmatrix <- asreml(fixed= flock.size ~ year,
                           random=~ vm(id, ped_inv) + ide(id) + vm(id2, natalenv_matrix),    ####
                           residual=~idv(units),
                           data=natalflocks,
                           na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                           )

summary(fsheritenvmatrix)$varcomp
vpredict(fsheritenvmatrix, envsimilarity ~ V1/(V1+V2+V3+V4)) #0.01928949 0.01365387
vpredict(fsheritenvmatrix, id ~ V3/(V1+V2+V3+V4)) #0.1825497 0.01385776
vpredict(fsheritenvmatrix, herit ~ V2/(V1+V2+V3+V4)) #0.01372072 0.01295723

#both matrices (Model i)

fsheritbothmatrix <- asreml(fixed= flock.size ~ year,
                            random=~ vm(id, ped_inv) + ide(id) + vm(id2, natalenv_matrix) + vm(id2, natalspat_matrix),    ####
                            residual=~idv(units),
                            data=natalflocks,
                            na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                            )
summary(fsheritbothmatrix)$varcomp 
vpredict(fsheritbothmatrix, repeatabilityid ~ V4/(V1+V2+V3+V4+V5)) #0.191272 0.01294572
vpredict(fsheritbothmatrix, herit ~ V3/(V1+V2+V3+V4+V5)) #0.0004429557 0.01083575
vpredict(fsheritbothmatrix, spatmatrix ~ V2/(V1+V2+V3+V4+V5)) #0.02717382 0.01611464
vpredict(fsheritbothmatrix, envsimmatrix ~ V1/(V1+V2+V3+V4+V5)) #9.416149e-09 1.726662e-10

1- pchisq(2 * (fsheritspatialmatrix$loglik - natalallheritonly$loglik),1) # 4.149565e-06 #adding the spatial matrix significantly improves the model
1- pchisq(2 * (fsheritenvmatrix$loglik - natalallheritonly$loglik),1) #0.001696816 #adding the natal env similarity matrix significantly improves the model

1- pchisq(2 * (fsheritbothmatrix$loglik - natalallheritonly$loglik),2) #2.49848e-05 #adding both the matrices significantly improves the model


#----------------------------------------------------------------------------------------------

#### ------- Degree natal models - all 3 winters with year as fixed effect ------

#only id (Model a)
degbaseid <- asreml(fixed= degree ~ year,  #### only id as random effect
                    random=~id,
                    residual=~idv(units),
                    data=natalwkndflocks,
                    na.action=na.method(x="omit", y="omit"))
summary(degbaseid)$varcomp
vpredict(degbaseid, repeatabilityid ~ V1/(V1+V2)) #0.415 (0.014)

#herit (Model b)
degallheritonly <- asreml(fixed= degree ~ year,
                          random=~ vm(id, ped_inv) + ide(id),    ####
                          residual=~idv(units),
                          data=natalwkndflocks,
                          na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                          maxit = 20)
summary(degallheritonly)$varcomp
vpredict(degallheritonly, herit ~ V1/(V1+V2+V3)) #0.041 (0.028)
vpredict(degallheritonly, id ~ (V1+V2)/(V1+V2+V3)) #0.416 (0.014)

#natal section (Model c)
degsectall <- asreml(fixed= degree ~ year,  #### id and section as re
                     random=~id + Section,
                     residual=~idv(units),
                     data=natalwkndflocks,
                     na.action=na.method(x="omit", y="omit"))

summary(degsectall)$varcomp
vpredict(degsectall, sect ~ V1/(V1+V2+V3)) #0.021 (0.013)
vpredict(degsectall, id ~ V2/(V1+V2+V3)) #0.388 (0.015)

#brood identity (Model d)
degpnumall <- asreml(fixed= degree ~ year,  #### id and pnum as re
                     random=~id + Pnum,
                     residual=~idv(units),
                     data=natalwkndflocks,
                     na.action=na.method(x="omit", y="omit"))
summary(degpnumall)$varcomp
vpredict(degpnumall, pnum ~ V1/(V1+V2+V3)) #0.043 (0.023)
vpredict(degpnumall, id ~ V2/(V1+V2+V3)) #0.372 (0.025)

#both natal effects (Model e)
degbothall <- asreml(fixed= degree ~ year,  #### id pnum and section as re
                     random=~id + Pnum + Section,
                     residual=~idv(units),
                     data=natalwkndflocks,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
                     maxit = 10)
summary(degbothall)$varcomp
vpredict(degbothall, sect ~ V1/(V1+V2+V3+V4)) #0.0205 (0.013)
vpredict(degbothall, pnum ~ V2/(V1+V2+V3+V4)) #0.0215 (0.022)
vpredict(degbothall, id ~ V3/(V1+V2+V3+V4)) #0.367 (0.025)

#both natal effects + pedigree (Model f)
degheritbothall <- asreml(fixed= degree ~ year,
                          random=~ vm(id, ped_inv) + ide(id) + Section + Pnum,    ####
                          residual=~idv(units),
                          data=natalwkndflocks,
                          na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                          maxit = 20)
summary(degheritbothall)$varcomp
vpredict(degheritbothall, sect ~ (V1)/(V1+V2+V3+V4+V5)) #0.0205 (0.013)
vpredict(degheritbothall, pnum ~ (V2)/(V1+V2+V3+V4+V5)) #0.0215 (0.022)
vpredict(degheritbothall, herit ~ (V3)/(V1+V2+V3+V4+V5)) #10^-8 (10^-10)
vpredict(degheritbothall, id ~ (V3+V4)/(V1+V2+V3+V4+V5)) #0.367 (0.025)

##Testing significance of random effects
degsectall$loglik #-8076.129
degpnumall$loglik #-8077.173
degbothall$loglik
degbaseid$loglik

#sect - id
1- pchisq(2 * (degsectall$loglik - degbaseid$loglik),1) #0 sect significant
#pnum - id
1- pchisq(2 * (degpnumall$loglik - degbaseid$loglik),1) #0 Pnum significant

#both - sect for pnum
1- pchisq(2 * (degbothall$loglik - degsectall$loglik),1) #0.148 Pnum not significant
1- pchisq(2 * (degbothall$loglik - degpnumall$loglik),1) #0 sect significant

1- pchisq(2 * (degbothall$loglik - degbaseid$loglik),2) #0 both effects significantly improve the model

#for models with matrices

# create another column to associate with second matrix

natalwkndflocks$id2 <- natalwkndflocks$id

#pedigree, id and spatial matrix (Model g)
degheritspatialmatrix <- asreml(fixed= degree ~ year,
                                random=~ vm(id, ped_inv) + ide(id) + vm(id2, natalspat_matrix),    ####
                                residual=~idv(units),
                                data=natalwkndflocks,
                                na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                                maxit = 20)

summary(degheritspatialmatrix)$varcomp
vpredict(degheritspatialmatrix, spatialdistance ~ V1/(V1+V2+V3+V4)) #0.07479738 0.03907895
vpredict(degheritspatialmatrix, id ~ V3/(V1+V2+V3+V4)) # 0.365566 0.021643
vpredict(degheritspatialmatrix, herit ~ V2/(V1+V2+V3+V4)) #5.308009e-08 2.421269e-09

#pedigree, id and natal env similarity matrix (Model h)
degheritenvmatrix <- asreml(fixed= degree ~ year,
                            random=~ vm(id, ped_inv) + ide(id) + vm(id2, natalenv_matrix),    ####
                            residual=~idv(units),
                            data=natalwkndflocks,
                            na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                            maxit = 20)

summary(degheritenvmatrix)$varcomp
vpredict(degheritenvmatrix, envsimilarity ~ V1/(V1+V2+V3+V4)) #0.06448233 0.0337805
vpredict(degheritenvmatrix, id ~ V3/(V1+V2+V3+V4)) # 0.3716263 0.02096789
vpredict(degheritenvmatrix, herit ~ V2/(V1+V2+V3+V4)) #1.38096e-07 5.410386e-09

#both matrices (Model i)
degheritbothmatrix <- asreml(fixed= degree ~ year,
                             random=~ vm(id, ped_inv) + ide(id) + vm(id2, natalenv_matrix) + vm(id2, natalspat_matrix),    ####
                             residual=~idv(units),
                             data=natalwkndflocks,
                             na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
                             maxit = 20)
summary(degheritbothmatrix)$varcomp 
vpredict(degheritbothmatrix, repeatabilityid ~ V4/(V1+V2+V3+V4+V5)) #0.365566 0.021643
vpredict(degheritbothmatrix, herit ~ V3/(V1+V2+V3+V4+V5)) #5.258138e-08 2.39852e-09
vpredict(degheritbothmatrix, spatmatrix ~ V2/(V1+V2+V3+V4+V5)) #0.07479734 0.03907895
vpredict(degheritbothmatrix, envsimmatrix ~ V1/(V1+V2+V3+V4+V5)) #3.769763e-08 1.719592e-09

1- pchisq(2 * (degheritspatialmatrix$loglik - degallheritonly$loglik),1) #8.437496e-10 adding the spatial matrix significantly improves the model
1- pchisq(2 * (degheritenvmatrix$loglik - degallheritonly$loglik),1) #1.418909e-07 adding the natal env similarity matrix significantly improves the model

1- pchisq(2 * (degheritbothmatrix$loglik - degallheritonly$loglik),2) #6.653421e-09 adding both the matrices significantly improves the model



# ---------- FIGURES ----------------

#Figure 5: Variance components from natal models for group size and degree

#Read in file with saved estimates from models 

datafsnatal <- read.csv("/ENTERYOURDIRECTORY/groupsizenatal_vcomp.csv", na.strings=c("", "NA"))
ggplot(datafsnatal, aes(x = model, y = Proportion, fill = Component)) +
  geom_bar(position = "stack",stat = "identity", show.legend=FALSE) +
  labs(x = "Models", y = "Proportion of variance", fill = "Variance Components") +
  ggtitle("Group Size")+
  scale_fill_manual(values = c(Vid = "#fabc60", Va = "#d96459", Vcsect = "black", Vbi = "#e276b5", Vlog = "#5FA683", Vr = "#40768B", Vcspat = "#58508d", Vcenv = "#60afc2"))+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")) + coord_cartesian(ylim = c(0.5, 1))+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), text=element_text(size=25), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
        plot.title = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0))) 


datadegnatal <- read.csv("/ENTERYOURDIRECTORY/degreenatal_vcomp.csv", na.strings=c("", "NA"))
ggplot(datadegnatal, aes(x = model, y = Proportion, fill = Component)) +
  geom_bar(position = "stack",stat = "identity", show.legend=FALSE) +
  labs(x = "Models", y = "Proportion of variance", fill = "Variance Components") +
  ggtitle("Degree")+
  scale_fill_manual(values = c(Vid = "#fabc60", Va = "#d96459", Vcsect = "black", Vbi = "#e276b5", Vr = "#40768B", Vcspat = "#58508d", Vcenv = "#60afc2"))+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")) + coord_cartesian(ylim = c(0.5, 1))+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), text=element_text(size=25), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
        plot.title = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0))) 


