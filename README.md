
This file contains descriptions of the data files provided with the code. Only data explicitly required for running all the models have been provided. These data come from the long-term individual-based study of great tits in Wytham Woods, Oxfordshire, UK. 


***Following are the focal datasets (1 and 2) containing individual records of social phenotypes. All the models use these phenotypes as response variables.***

(1) rawflockingevents.rds and flockingevents_natalinfo.rds 
__________________________________________________________

id: unique identifier for each individual bird observed at feeders
sex: sex of individual
adjuv: age class of individual
born: whether locally-born (TRUE) or immigrant (FALSE)
logger: unique identifier for RFID logger at feeders
year: the year (winter) in which flocking event (observation) occured (2011-2013)
nwkend: the weekend number during the winter when the flocking event occured (1-13)
date: date of flocking event
time: time of flocking event
flock: unique identifier for the flocking event
flock.size: number of individuals observed in a flocking event (group size) 

Section: natal section, region where individual was born

(Lay.date, April.lay.date, Hatch.date, April.hatch.date, Incubation.duration, Clutch.size, Num.chicks, April.ring.date, Mean.chick.weight - variables not used for any of the analyses reported here)

nestbox: the nestbox in which individual was born
Pnum: the brood identity or a unique identifier for the brood an individual belongs to


(2) socialnetworkdatafull.rds and socialnetwork_natalinfo.rds
_____________________________________________________________

id: unique identifier for each individual bird observed at feeders
max.flock.size: maximum groupsize of the individual during the weekend (sampling period)
mean.flock.size: average groupsize of the individual over all the flocking events in the weekend (sampling period)
nr.flocks: number of flocks the individual was part of throughout sampling period
nwkend: the weekend number (social networks generated for each weekend in every winter), and all metrics thus correspond to a particular weekend 
degree: number of associations in the network
w.degree: strength, number of associations in the network weighted by frequency
eigen_cent: eigenvector centrality of individual in the network
w.eigen_cent: weighted eigenvector centrality of individual in the network
betweenness: betweenness of individual in the network
year: the winter year in which network created for each weekend in it
sex: sex of individual
adjuv: age-class of individual
born: whether locally-born (TRUE) or immigrant (FALSE)

Natal variables in (2) is the same as in (1)


(3) GTITprunedpedigree.csv
___________________________

the pedigree is a pruned version of the overall pedigree from 1960-2022, containing only informative individuals for the analyses

id: unique identifier for each bird in the system
dam: mother id
sire: father id

(4) natalhabitatdata.csv
_________________________

pnum: the brood identity or a unique identifier for the brood an individual belongs to
id, Section, nestbox, birthyear same as before

x and y: coordinates of natal nestbox
edge_edi: distance of nestbox from the edge of the woodland
altitude_m: altitude of natal nestbox
northness: northness of the natal nestbox
no_trees_75m: number of oak trees within 75 metres of natal nestbox
area_polygon: area of the territory polygon around natal nestbox
area_polygon_sqrt: square root of area_polygon

***Following datasets contain the outputs (variance componenets, estimates and standard error) from the models. These have been used to make the plots.***

(5) groupsize model estimate.csv and socnet model estimates.csv
________________________________________________________________

Year: winter year
Model: Model number (detailed in methods)
EstimateType: whether the estimate is for repeatability or heritability
Estimate: estimate obtained from model output
SE: standar error of estimate
Trait: social network trait

(6) groupsizenatal_vcomp.csv and degreenatal_vcomp.csv
_______________________________________________________

model: model number (detailed in methods)
Vp: total sum of variance components (total phenotypic variance)
Component: type of variance component
Variance: Variance component (num)
Proportion: Amount of variance out of total phenotypic variance
Percentage: Percentage of variance out of total phenotypic variance


Code files:

Groupsize animal models.R
Social Network Traits animal models.R
Natal Effects animal models.R


