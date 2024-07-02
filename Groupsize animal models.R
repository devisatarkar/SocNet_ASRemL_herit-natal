
#### ANIMAL MODELS FOR GROUP SIZE #####

#REQUIRED LIBRARIES

library(dplyr)
library(asreml)
library(ggplot2)

# ---------- PREPARE DATA ----------------

#Read in source data (all flocking events)

ogflocks <- readRDS("/ENTERYOURDIRECTORY/rawflockingevents.rds") 

#remove experimental loggers
values_to_exclude <- c("5A1", "5A2", "5B1", "5B2", "5C1", "5C2", "8C1", "8C2", "8D1", "8D2", "8E1", "8E2")

flocks <- ogflocks %>%
  filter(!logger %in% values_to_exclude)
print(length(unique(flocks$id))) #There are 1823 unique individuals (791818 flocking events)

#Read in pedigree data 

ped_data <- read.csv("/ENTERYOURDIRECTORY/GTITprunedpedigree.csv", na.strings=c("", "NA"))

#Make data usable for ASReml

flocks <- flocks %>% dplyr::select(id, sex, adjuv, born, logger, year, nwkend, flock, flock.size)
flocks <- flocks %>% dplyr::mutate_at(c('id', 'sex', 'adjuv', 'born', 'logger', 'year', 'nwkend', 'flock'), as.factor) 
flocks$id <- factor(toupper(levels(flocks$id))[flocks$id])

ped_inv <- ainverse(ped_data) #inverse of pedigree data - needed for ASReml

# ---------- RUN MODELS ----------------

#2011

flocks2011 <- flocks %>% filter(year == 2011)
print(length(unique(flocks2011$id))) #1085 (343458 obs)
print(length(unique(flocks2011$logger))) #65

#--------------------------------------------------------------------------------
model12011 <- asreml(fixed= flock.size ~ 1,  #### only id as random effect
                     random=~id,
                     residual=~idv(units),
                     data=flocks2011,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model12011)$varcomp
vpredict(model12011, repeatabilityid ~ V1/(V1+V2)) #0.192 (0.007)


model22011 <- asreml(fixed= flock.size ~ 1,                   
                     random=~ vm(id, ped_inv) + ide(id),    #### ped + id random effect
                     residual=~idv(units),
                     data=flocks2011,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model22011)$varcomp
vpredict(model22011, repeatabilityid ~ (V1+V2)/(V1+V2+V3)) #0.192 (0.007)
vpredict(model22011, h2 ~ (V1)/(V1+V2+V3)) #2.523159e-08 

1 - pchisq(2 * (model22011$loglik - model12011$loglik), 1) #1 #adding the pedigree random effect does not improve the model

model32011 <- asreml(fixed= flock.size ~ 1,
                     random=~logger,            ####logger as random effect
                     residual=~idv(units),
                     data=flocks2011,
                     na.action=na.method(x="omit", y="omit"),workspace = 64e6,  # Boost memory
)
summary(model32011)$varcomp
vpredict(model32011, repeatabilitylogger ~ (V1)/(V1+V2)) #0.292 (0.036)


model42011 <- asreml(fixed= flock.size ~ 1,
                     random=~id + logger,          ####id and logger as random effect
                     residual=~idv(units),
                     data=flocks2011,
                     na.action=na.method(x="omit", y="omit"), workspace = 32e6,  # Boost memory
)
summary(model42011)$varcomp
vpredict(model42011, repeatabilitylogger ~ (V1)/(V1+V2+V3)) #0.318 (0.038)
vpredict(model42011, repeatabilityid ~ (V2)/(V1+V2+V3)) #0.048 (0.003)


1 - pchisq(2 * (model42011$loglik - model12011$loglik), 1) #0 #adding the logger term significantly improves the model

model52011 <- asreml(fixed= flock.size ~ 1,
                     random=~ vm(id, ped_inv) + ide(id) + logger,    ####id, logger and ped as random effects
                     residual=~idv(units),
                     data=flocks2011,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model52011)$varcomp
vpredict(model52011, repeatabilitylogger ~ (V1)/(V1+V2+V3+V4)) #0.317 (0.038)
vpredict(model52011, repeatabilityid ~ (V2+V3)/(V1+V2+V3+V4)) #0.05 (0.003)
vpredict(model52011, h2 ~ (V2)/(V1+V2+V3+V4)) #0.004 (0.003)

#2012

flocks2012 <- flocks %>% filter(year == 2012)
print(length(unique(flocks2012$id))) #720 (244064 obs)
print(length(unique(flocks2012$logger))) #65

#--------------------------------------------------------------------------------

model12012 <- asreml(fixed= flock.size ~ 1,  #### only id as random effect
                     random=~id,
                     residual=~idv(units),
                     data=flocks2012,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model12012)$varcomp
vpredict(model12012, repeatabilityid ~ V1/(V1+V2)) #0.303 (0.011)


model22012 <- asreml(fixed= flock.size ~ 1,                   
                     random=~ vm(id, ped_inv) + ide(id),    #### ped + id random effect
                     residual=~idv(units),
                     data=flocks2012,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model22012)$varcomp
vpredict(model22012, repeatabilityid ~ (V1+V2)/(V1+V2+V3)) #0.306
vpredict(model22012, h2 ~ (V1)/(V1+V2+V3)) #0.019

model32012 <- asreml(fixed= flock.size ~ 1,
                     random=~logger,
                     residual=~idv(units),
                     data=flocks2012,
                     na.action=na.method(x="omit", y="omit"),workspace = 64e6,  # Boost memory
)
summary(model32012)$varcomp
vpredict(model32012, repeatabilitylogger ~ (V1)/(V1+V2)) #0.296


model42012 <- asreml(fixed= flock.size ~ 1,
                     random=~id + logger,
                     residual=~idv(units),
                     data=flocks2012,
                     na.action=na.method(x="omit", y="omit"), workspace = 32e6,  # Boost memory
)
summary(model42012)$varcomp
vpredict(model42012, repeatabilitylogger ~ (V1)/(V1+V2+V3)) #0.332
vpredict(model42012, repeatabilityid ~ (V2)/(V1+V2+V3)) #0.045


model52012 <- asreml(fixed= flock.size ~ 1,
                     random=~ vm(id, ped_inv) + ide(id) + logger,    ####
                     residual=~idv(units),
                     data=flocks2012,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model52012)$varcomp
vpredict(model52012, repeatabilitylogger ~ (V1)/(V1+V2+V3+V4)) #0.332
vpredict(model52012, repeatabilityid ~ (V2+V3)/(V1+V2+V3+V4)) #0.045
vpredict(model52012, h2 ~ (V2)/(V1+V2+V3+V4)) #10^-8


#2013

flocks2013 <- flocks %>% filter(year == 2013)
print(length(unique(flocks2013$id))) #789 (204296 obs)
print(length(unique(flocks2013$logger))) #65

#--------------------------------------------------------------------------------

model12013 <- asreml(fixed= flock.size ~ 1,  #### only id as random effect
                     random=~id,
                     residual=~idv(units),
                     data=flocks2013,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model12013)$varcomp
vpredict(model12013, repeatabilityid ~ V1/(V1+V2)) #0.205 (0.008)


model22013 <- asreml(fixed= flock.size ~ 1,                   
                     random=~ vm(id, ped_inv) + ide(id),    #### ped + id random effect
                     residual=~idv(units),
                     data=flocks2013,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model22013)$varcomp
vpredict(model22013, repeatabilityid ~ (V1+V2)/(V1+V2+V3)) #0.209 (0.009)
vpredict(model22013, h2 ~ (V1)/(V1+V2+V3)) #0.025 (0.015)


model32013 <- asreml(fixed= flock.size ~ 1,
                     random=~logger,
                     residual=~idv(units),
                     data=flocks2013,
                     na.action=na.method(x="omit", y="omit"),workspace = 64e6,  # Boost memory
)
summary(model32013)$varcomp
vpredict(model32013, repeatabilitylogger ~ (V1)/(V1+V2)) #0.279


model42013 <- asreml(fixed= flock.size ~ 1,
                     random=~id + logger,
                     residual=~idv(units),
                     data=flocks2013,
                     na.action=na.method(x="omit", y="omit"), workspace = 32e6,  # Boost memory
)
summary(model42013)$varcomp
vpredict(model42013, repeatabilitylogger ~ (V1)/(V1+V2+V3)) #0.294
vpredict(model42013, repeatabilityid ~ (V2)/(V1+V2+V3)) #0.03


model52013 <- asreml(fixed= flock.size ~ 1,
                     random=~ vm(id, ped_inv) + ide(id) + logger,    ####
                     residual=~idv(units),
                     data=flocks2013,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model52013)$varcomp
vpredict(model52013, repeatabilitylogger ~ (V1)/(V1+V2+V3+V4)) #0.294
vpredict(model52013, repeatabilityid ~ (V2+V3)/(V1+V2+V3+V4)) #0.030
vpredict(model52013, h2 ~ (V2)/(V1+V2+V3+V4)) #10^-8


#Between-year (all 3 winters)

#year as fixed effect
model1all <- asreml(fixed= flock.size ~ year,  #### only id as random effect
                    random=~id,
                    residual=~idv(units),
                    data=flocks,
                    na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model1all)$varcomp
vpredict(model1all, repeatabilityid ~ V1/(V1+V2)) #0.209 (0.005)
summary(model1all, coef = TRUE)$coef.fixed
wald.asreml(model1all, ssType = "conditional", denDF = "numeric")


#year as random effect *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_

model1allre <- asreml(fixed= flock.size ~ 1,  #### only id as random effect
                      random=~id + year,
                      residual=~idv(units),
                      data=flocks,
                      na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model1allre)$varcomp
vpredict(model1allre, repeatabilityid ~ V2/(V1+V2+V3)) #0.199 (0.0099)

#base model for random effect *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*

model1all0 <- asreml(fixed= flock.size ~ 1,  #### only id as random effect
                     random=~id,
                     residual=~idv(units),
                     data=flocks,
                     na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model1all0)$varcomp
vpredict(model1all0, repeatabilityid ~ V1/(V1+V2)) #0.230 (0.006)

1 - pchisq(2 * (model1allre$loglik - model1all0$loglik), 1) #0 


#using year as fixed effect in all future models ---------------------------

model2all <- asreml(fixed= flock.size ~ year,                   
                    random=~ vm(id, ped_inv) + ide(id),    #### ped + id random effect
                    residual=~idv(units),
                    data=flocks,
                    na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
)
summary(model2all)$varcomp
vpredict(model2all, repeatabilityid ~ (V1+V2)/(V1+V2+V3)) #0.208 (0.007)
vpredict(model2all, h2 ~ (V1)/(V1+V2+V3)) #0.002 (0.007)


model3all <- asreml(fixed= flock.size ~ year,
                    random=~logger,
                    residual=~idv(units),
                    data=flocks,
                    na.action=na.method(x="omit", y="omit"),workspace = 64e6,  # Boost memory
)
summary(model3all)$varcomp
vpredict(model3all, repeatabilitylogger ~ (V1)/(V1+V2)) #0.245


model4all <- asreml(fixed= flock.size ~ year,
                    random=~id + logger,
                    residual=~idv(units),
                    data=flocks,
                    na.action=na.method(x="omit", y="omit"), workspace = 64e6,  # Boost memory
)
summary(model4all)$varcomp
vpredict(model4all, repeatabilitylogger ~ (V1)/(V1+V2+V3)) #0.261
vpredict(model4all, repeatabilityid ~ (V2)/(V1+V2+V3)) #0.082


model5all <- asreml(fixed= flock.size ~ year,
                    random=~ vm(id, ped_inv) + ide(id) + logger,    ####
                    residual=~idv(units),
                    data=flocks,
                    na.action=na.method(x="omit", y="omit"), workspace = 96e6,  # Boost memory
)
summary(model5all)$varcomp
vpredict(model5all, repeatabilitylogger ~ (V1)/(V1+V2+V3+V4)) #0.261
vpredict(model5all, repeatabilityid ~ (V2+V3)/(V1+V2+V3+V4)) #0.08
vpredict(model5all, h2 ~ (V2)/(V1+V2+V3+V4)) #0.001

summary(model5all, coef = TRUE)$coef.fixed
wald.asreml(model5all, ssType = "conditional", denDF = "numeric")

# ---------- FIGURES ----------------

#FIGURE 2: Repeatability and Heritability estimates from Groupsize models

#Read in file with saved estimates from models 
fsdata <- read.csv("/ENTERYOURDIRECTORY/groupsize model estimates.csv", na.strings=c("", "NA")) 
fsheritdata <- fsdata %>% filter(Model=="Model 2")
fsloggerdata <- fsdata %>% filter(!Model=="Model 2")


fsplot1 <-ggplot(fsheritdata, aes(y = Year, x = Estimate, color = EstimateType)) +
  geom_point(size = 4, shape = 16, aes(fill = EstimateType)) +
  geom_errorbarh(aes(xmin = Estimate - SE, xmax = Estimate + SE),
                 height = 0.1, color="black" ) +
  labs(y = "Year", x = "Estimate") +
  ggtitle("Repeatability & Heritability of Group Size") +
  scale_fill_manual(values = c("#d96459", "#1FBBC6")) +
  scale_shape_manual(values = c(16, 17)) + 
  theme(legend.position = "bottom",  # Change legend position
        legend.title = element_text(size = 16),  # Change legend title font size
        legend.text = element_text(size = 14),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))+
  scale_color_manual(values = c("#d96459", "#1FBBC6")) + 
  theme(axis.text=element_text(size=20), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        plot.title = element_text(size=20, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        panel.grid.major.y = element_line(color = "grey", linetype = "dotted")) 

fsplot2 <- ggplot(fsloggerdata, aes(y = Year, x = Estimate, color = EstimateType, group = Model)) +
  geom_point(position = position_dodge(width = 0.9), size = 4) +
  geom_errorbarh(aes(xmin = Estimate - SE, xmax = Estimate + SE),
                 position = position_dodge(width = 0.9), width = 0.05, color="black", height=0.1) +
  labs(y = "Year", x = "Estimate", color = "Estimate Type", shape = "Model") +
  ggtitle("Repeatability of Group Size with spatial effects") +
  scale_fill_manual(values = c("#316B4E", "#1FBBC6")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_y_discrete(labels = fsloggerdata$Year,  # Use Year as labels
                   breaks = fsloggerdata$Year) +  # Set breaks to match Year values
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove major grid lines on x-axis
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("#1FBBC6", "#316B4E")) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title.x = element_text(size = 20, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        plot.title = element_text(size = 20, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        panel.grid.major.y = element_line(color = "grey", linetype = "dotted"))  # Add horizontal grid lines


#FIGURE 3: Variance Components 

#Following is the plot for 2011 models. Replace the 2011 models with the models from other years to obtain plots for the corresponding winters

id2011 <- c(summary(model12011)$varcomp[1,"component"],0,0, summary(model12011)$varcomp[2,"component"],
            summary(model12011)$varcomp[1,"component"]+summary(model12011)$varcomp[2,"component"])

idped2011 <- c(summary(model22011)$varcomp[2,"component"],
               summary(model22011)$varcomp[1,"component"],0,summary(model22011)$varcomp[3,"component"],
               summary(model22011)$varcomp[1,"component"]+ summary(model22011)$varcomp[2,"component"]+summary(model22011)$varcomp[3,"component"])

logger2011 <- c(0,0,summary(model32011)$varcomp[1,"component"],summary(model32011)$varcomp[2,"component"],
                summary(model32011)$varcomp[1,"component"]+summary(model32011)$varcomp[2,"component"])

idlog2011 <- c(summary(model42011)$varcomp[2,"component"],0,
               summary(model42011)$varcomp[1,"component"],summary(model42011)$varcomp[3,"component"],
               summary(model42011)$varcomp[1,"component"]+ summary(model42011)$varcomp[2,"component"]+summary(model42011)$varcomp[3,"component"])

idlogped2011 <- c(summary(model52011)$varcomp[3,"component"],summary(model52011)$varcomp[2,"component"],
                  summary(model52011)$varcomp[1,"component"], summary(model52011)$varcomp[4,"component"],
                  summary(model52011)$varcomp[1,"component"]+
                    summary(model52011)$varcomp[2,"component"]+
                    summary(model52011)$varcomp[3,"component"]+
                    summary(model52011)$varcomp[4,"component"])

id2011data <- data.frame(id2011, idped2011,logger2011,idlog2011, idlogped2011)
rownames(id2011data) <- c("Vid", "Va","Vlog","Vr", "Vp")

id2011data <- t(id2011data)

id2011data <- cbind(model = model_values, id2011data)
id2011df <- as.data.frame(id2011data)
id2011df <- id2011df %>%
  mutate(across(starts_with("V"), as.numeric))

data2011long <- id2011df %>%
  pivot_longer(cols = c(Vid, Va, Vlog, Vr), names_to = "Component", values_to = "Variance")
data2011long <- data2011long %>%
  group_by(model) %>%
  mutate(Proportion = Variance / sum(Variance),
         Percentage = Proportion * 100)

barplot2011 <- ggplot(data2011long, aes(x = model, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", show.legend=FALSE) +
  labs(x = "Models", y = "Proportion of variance", fill = "Variance Components") +
  ggtitle("2011 winter")+
  scale_fill_manual(values = c(Vid = "#fabc60", Va = "#d96459", Vlog = "#5FA683", Vr = "#40768B"))+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"))


