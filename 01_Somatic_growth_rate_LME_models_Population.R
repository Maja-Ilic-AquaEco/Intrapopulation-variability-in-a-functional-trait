############################################################################################
##                                                                                        ##
## R code from the manuscript "Intrapopulation variability in a functional trait:         ##
## susceptibility of Daphnia to limitation by dietary fatty acids"                        ##
##                                                                                        ##
## Authors: Maja Ilic (1, 2, a), Mathilde Cordellier (3, b) and Patrick Fink (1, 4, 5, c) ##  
##                                                                                        ##      
## (1) University of Cologne, Institute for Zoology,                                      ##  
##     Zülpicher Strasse 47b, 50674 Köln, Germany                                         ##
## (2) Present address: Queen's University Belfast, School of Biological Sciences,        ##
##     19 Chlorine Gardens, Belfast BT9 5DL, United Kingdom                               ##
## (3) Universität Hamburg, Institute of Zoology,                                         ##
##     Martin-Luther-King Platz 3, 20146 Hamburg, Germany                                 ##
## (4) Present address: Helmholtz Centre of Environmental Research - UFZ,                 ##
##     Department River Ecology, Brückstrasse 3a, 39114 Magdeburg, Germany AND            ##
## (5) Helmholtz Centre of Environmental Research - UFZ, Department Aquatic Ecosystem     ##
##     Analysis, Brückstrasse 3a, 39114 Magdeburg, Germany                                ##
##                                                                                        ##
## E-Mail:                                                                                ##
## (a) Author for correspondance: maja.ilic.bio@gmail.com,                                ##
##                                ORCID iD: 0000-0002-8387-9932                           ##
## (b) mathilde.cordellier@uni-hamburg.de                                                 ##
## (c) patrick.fink@ufz.de                                                                ##
##                                                                                        ##
##                                                                                        ##
## Part 1: Somatic growth rate of Daphnia as response variable (entire population)        ##  
##                                                                                        ##
## Containts code for data import and statistical analyses:                               ##
## Linear mixed-effects model - development and validation following:                     ##
##                                                                                        ##
## Zuur, A.F., Ieno, E.N., Walker, N.J., Saveliev, A.A., Smith, G.M., 2009.               ##
## Mixed effects models and extensions in ecology with R, Statistics for biology          ##
## and health. Springer, New York, NY.                                                    ##
##                                                                                        ##    
############################################################################################

##===============================================
## Packages ----

library(readxl)    ## for data import from an Excel-file, read.excel()
library(car)       ## for leveneTest() and Anova()
library(nlme)      ## for linear mixed-effects models, lme()
library(lsmeans)   ## for pairwise comparisons of the groups from the model lsmeans()
library(ape)       ## for an easier access of variance within the model, varcomp()
library(dplyr)     ## for data summary

##===============================================
## Data import ----

myfile <- read_excel("Raw data Ilic, Cordellier and Fink, Freshwater Biology.xlsx", sheet = "SGR")
head(myfile, 10)  ## first 10 lines of the imported dataset
str(myfile)       ## structure of the dataset and classes of the variables within the dataset

########################################################################################
##                                                                                    ##    
## Explanation of the variables:                                                      ##
##                                                                                    ##
## Day: Gives a day of the year 2017 on which each individual experiment was started  ##
##      within the experimental period (March - November 2017)                        ##
##      Day 1: 1st of January, 2017                                                   ##
##                                                                                    ##
## Experiment_Nr: Number of the experiment for every of the genotypes                 ##
##                (at least 3 experiments per genotype)                               ##
##                                                                                    ##
## Genotype: Laboratory code/name for each D. longispina genotype isolated from       ##
##           the Lake Klostersee (KL) and used in this study                          ##
##                                                                                    ##
## Species: Species, identified following Rusek et al. (2015).                        ##
##          Note that 11 out of 12 genotypes are D. longispina, while only            ##
##          the genotype KL8 is a D. longispina × cucullata hybrid                    ##
##                                                                                    ##
## PUFA: Liposomes used to supplement the diet (green alga Acutodesmus obliquus)      ##
##       with following PUFAs: C - control liposomes, PUFA-free;                      ##
##                             ALA - alpha-linolenic acid;                            ##
##                             EPA - eicosapentaenoic acid;                           ##
##                             ARA - arachidonic acid                                 ##
##                                                                                    ##
## Growth_rate: Somatic growth rate of Daphnia individuals                            ##
##              (after 6 days of experiment) in every jar per treatment               ##
##              (see Liposome), calculated with the Eq. 1 (see M and M)               ##
##                                                                                    ##
########################################################################################

##===============================================
## Specify a certain order of the factors PUFA and Genotype (not mandatory) ----

myfile$PUFA <- factor(myfile$PUFA, levels = c("C","ALA","EPA","ARA"))

myfile$Genotype <- factor(myfile$Genotype, levels = c("KL14","KL3","KL8","KL11",
                                                      "KL13","KL50","KL53","KL54",
                                                      "KL73","KL82","KL83","KL93"))

nr_exp <- with(myfile,tapply(Experiment_Nr,list(Genotype),max))
nr_exp ## gives the number of experiments performed per genotype

#######################################
##                                   ##
##     All 12 genotypes included     ##
##                                   ##
#######################################

##===============================================
## Check normality of the raw data ----

shapiro.test(myfile$Growth_rate)   ## normal distribution not given
par(mfrow = c(1,2))
hist(myfile$Growth_rate, las = 1)
qqPlot(myfile$Growth_rate, distribution = "norm", las = 1)

##===============================================
## Check the homogeneity of variances ----

leveneTest(Growth_rate ~ PUFA, data = myfile)          ## heterogeneious (unequal) variances
leveneTest(Growth_rate ~ Genotype, data = myfile)      ## heterogeneious variances
leveneTest(Growth_rate ~ PUFA*Genotype, data = myfile) ## heterogeneious variances

dev.off()
split.screen(c(2,1))
split.screen(c(1,2), screen = 1)
screen(2)
par(mar = c(6,5,1,1))
boxplot(Growth_rate ~ PUFA*Genotype, data = myfile, las = 2, xlab = "")
screen(3)
par(mar = c(4,5,1,1))
boxplot(Growth_rate ~ PUFA, data = myfile, las = 1)
screen(4)
par(mar = c(4,5,1,1))
boxplot(Growth_rate ~ Genotype, data = myfile, las = 2)

##===============================================
## Null model: linear mixed-effects model (LME) ----

# Response variable: Growth_rate (somatic growth rate of Daphnia)
# Fixed effects: PUFA
#                Genotype
#                Interaction PUFA x Genotype
# Random effect: Time, given in days (Day)
# Variance structure (Weights): not specified (NULL) 

Model.null <- lme(Growth_rate ~ PUFA*Genotype, 
                  random = ~1|Day,
                  data = myfile,
                  method = "REML")
(summ.null <- summary(Model.null))

##===============================================
## Model.1 ----
# Variance structure: different variances allowed per stratum (only PUFA)

vf1 <- varIdent(form= ~ 1 | PUFA)

Model.1 <- lme(Growth_rate ~ PUFA*Genotype, 
               random = ~1|Day,
               weights = vf1,
               data = myfile, 
               method = "REML")
(summ.1 <- summary(Model.1))

##===============================================
## Model.2 ----
# Variance structure: different variances allowed per stratum (only Genotype)

vf2 <- varIdent(form= ~ 1 | Genotype)

Model.2 <- lme(Growth_rate ~ PUFA*Genotype, 
               random = ~1|Day,
               weights = vf2,
               data = myfile,
               method = "REML")
(summ.2 <- summary(Model.2))

##===============================================
## Model.3
# Variance structure: different variances allowed per stratum (PUFA and Genotype)

vf3 <- varIdent(form= ~ 1 | PUFA*Genotype)

Model.3 <- lme(Growth_rate ~ PUFA*Genotype, 
               random = ~1|Day,
               weights = vf3,
               data = myfile, 
               method = "REML")
(summ.3 <- summary(Model.3))

##===============================================
## Model selection decided using Aikaike Information Criterion (AIC) ----
# The model with the lowest AIC was chosen for further statistical analysis

AIC(Model.null, Model.1, Model.2, Model.3)  ## lowest AIC: Model.3

##===============================================
## Anova (F and Chisq statistics) ----

anova(Model.3) 

##===============================================
## Model validation ----

E.3 <- resid(Model.3, type = "normalized")  ## residuals
F.3 <- fitted(Model.3)                      ## fitted values

shapiro.test(E.3) ## not normally distributed
dev.off()
par(mfrow = c(1,1))
hist(E.3,las = 1)

# Note: due to a large sample size (n = 489), the deviation from the normal (Gaussian) 
# distribution can be ignored (Underwood, 1997)

leveneTest(E.3 ~ PUFA, data = myfile)             ## homogeneous variances
leveneTest(E.3 ~ Genotype, data = myfile)         ## homogeneous variances
leveneTest(E.3 ~ PUFA*Genotype, data = myfile)    ## homogeneous variances

##===============================================
## Graphical Validation ----

min(E.3)
max(E.3)

dev.off()
split.screen(c(2,1))
split.screen(c(1,2), screen = 1)
screen(2)
par(mar = c(4,4,2,2))
boxplot(E.3 ~ Genotype, data = myfile, main = "Model.3", 
        xlab = "Genotype", ylab = "Residuals", las = 2, ylim = c(-3,3))
screen(3)
par(mar = c(4,4,2,2))
plot(x = F.3, y = E.3, main = "Model.3",
     xlab = "Fitted values", ylab = "Residuals", las = 1, ylim = c(-3,3))
screen(4)
par(mar = c(4,4,2,2))
boxplot(E.3 ~ PUFA, data = myfile, main = "Model.3", 
        xlab = "PUFA", ylab = "Residuals", las = 1, ylim = c(-3,3))

par(mfrow = c(1,2))
qqPlot(resid(Model.3, type = "normalized"), distribution = "norm", las = 1)
qqPlot(resid(Model.3, type = "response"), distribution = "norm", las = 1)

##===============================================
## Goodness-of-fit ----

cor(myfile$Growth_rate, fitted(Model.3))^2

##===============================================
## Pairwise comparisons between the groups ----
# Requires library(lsmeans)

options(max.print = 999999999)  ## to be able to display the output in its full length

lsmeans(Model.3, list(pairwise ~ PUFA), adjust = "tukey")
lsmeans(Model.3, list(pairwise ~ Genotype), adjust = "tukey")
lsmeans(Model.3, list(pairwise ~ PUFA*Genotype), adjust = "tukey")

##===============================================
## How much of the variance is explained by the fixed and random factors? ----
# Requires library(ape)

(v1 <- varcomp(Model.3, TRUE, FALSE)*100)  ## in percent (sums up to 100%)
(v2 <- varcomp(Model.3, FALSE, FALSE))

#####################################################
##                                                 ##
##      Only D. longispina genotypes included      ##
##  Hybrid KL8 excluded from statistical analysis  ##
##                                                 ##
#####################################################

##===============================================
## Select all genotypes that are D. longispina ----

mylongi <- subset(myfile, Species == "D. longispina")
mylongi <- droplevels(mylongi)

##===============================================
## Check normality of the raw data ----

shapiro.test(mylongi$Growth_rate)   ## normal distribution not given
dev.off()
par(mfrow = c(1,2))
hist(mylongi$Growth_rate, las = 1)
qqPlot(mylongi$Growth_rate, distribution = "norm", las = 1)

##===============================================
## Check the homogeneity of variances ----

leveneTest(Growth_rate ~ PUFA, data = mylongi)          ## heterogeneious (unequal) variances
leveneTest(Growth_rate ~ Genotype, data = mylongi)      ## heterogeneious variances
leveneTest(Growth_rate ~ PUFA*Genotype, data = mylongi) ## heterogeneious variances

dev.off()
split.screen(c(2,1))
split.screen(c(1,2), screen = 1)
screen(2)
par(mar = c(6,5,1,1))
boxplot(Growth_rate ~ PUFA*Genotype, data = mylongi, las = 2, xlab = "")
screen(3)
par(mar = c(4,5,1,1))
boxplot(Growth_rate ~ PUFA, data = mylongi, las = 1)
screen(4)
par(mar = c(4,5,1,1))
boxplot(Growth_rate ~ Genotype, data = mylongi, las = 2)

##===============================================
## Null model: linear mixed-effects model (LME) ----

# Response variable: Growth_rate (somatic growth rate of Daphnia)
# Fixed effects: PUFA
#                Genotype
#                Interaction PUFA x Genotype
# Random effect: Time, given in days (Day)
# Variance structure (Weights): not specified (NULL) 

Model.null <- lme(Growth_rate ~ PUFA*Genotype, 
                  random = ~1|Day,
                  data = mylongi,
                  method = "REML")
(summ.null <- summary(Model.null))

##===============================================
## Model.1 ----
# Variance structure: different variances allowed per stratum (only PUFA)

vf1 <- varIdent(form= ~ 1 | PUFA)

Model.1 <- lme(Growth_rate ~ PUFA*Genotype, 
               random = ~1|Day,
               weights = vf1,
               data = mylongi, 
               method = "REML")
(summ.1 <- summary(Model.1))

##===============================================
## Model.2 ----
# Variance structure: different variances allowed per stratum (only Genotype)

vf2 <- varIdent(form= ~ 1 | Genotype)

Model.2 <- lme(Growth_rate ~ PUFA*Genotype, 
               random = ~1|Day,
               weights = vf2,
               data = mylongi,
               method = "REML")
(summ.2 <- summary(Model.2))

##===============================================
## Model.3 ----
# Variance structure: different variances allowed per stratum (PUFA and Genotype)

vf3 <- varIdent(form= ~ 1 | PUFA*Genotype)

Model.3 <- lme(Growth_rate ~ PUFA*Genotype, 
               random = ~1|Day,
               weights = vf3,
               data = mylongi, 
               method = "REML")
(summ.3 <- summary(Model.3))

##===============================================
## Model selection decided using Aikaike Information Criterion (AIC) ----
# The model with the lowest AIC was chosen for further statistical analysis

AIC(Model.null, Model.1, Model.2, Model.3)  ## lowest AIC: Model.3

##===============================================
## Anova (F and Chisq statistics) ----

anova(Model.3) 

##===============================================
## Model validation ----

E.3 <- resid(Model.3, type = "normalized")  ## residuals
F.3 <- fitted(Model.3)                      ## fitted values

shapiro.test(E.3) ## not normally distributed
dev.off()
par(mfrow = c(1,1))
hist(E.3,las = 1)

# Note: due to a large sample size (n = 453), the deviation from the normal (Gaussian) 
# distribution can be ignored (Underwood, 1997)

leveneTest(E.3 ~ PUFA, data = mylongi)             ## homogeneous variances
leveneTest(E.3 ~ Genotype, data = mylongi)         ## homogeneous variances
leveneTest(E.3 ~ PUFA*Genotype, data = mylongi)    ## homogeneous variances

##===============================================
## Graphical Validation ----

min(E.3)
max(E.3)

dev.off()
split.screen(c(2,1))
split.screen(c(1,2), screen = 1)
screen(2)
par(mar = c(4,4,2,2))
boxplot(E.3 ~ Genotype, data = mylongi, main = "Model.3", 
        xlab = "Genotype", ylab = "Residuals", las = 2, ylim = c(-3,3))
screen(3)
par(mar = c(4,4,2,2))
plot(x = F.3, y = E.3, main = "Model.3",
     xlab = "Fitted values", ylab = "Residuals", las = 1, ylim = c(-3,3))
screen(4)
par(mar = c(4,4,2,2))
boxplot(E.3 ~ PUFA, data = mylongi, main = "Model.3", 
        xlab = "PUFA", ylab = "Residuals", las = 1, ylim = c(-3,3))

par(mfrow = c(1,2))
qqPlot(resid(Model.3, type = "normalized"), distribution = "norm", las = 1)
qqPlot(resid(Model.3, type = "response"), distribution = "norm", las = 1)

##===============================================
## Goodness-of-fit ----

cor(mylongi$Growth_rate, fitted(Model.3))^2

##===============================================
## Pairwise comparisons between the groups ----
# Requires library(lsmeans)

options(max.print = 999999999)  ## to be able to display the output in its full length

lsmeans(Model.3, list(pairwise ~ PUFA), adjust = "tukey")
lsmeans(Model.3, list(pairwise ~ Genotype), adjust = "tukey")
lsmeans(Model.3, list(pairwise ~ PUFA*Genotype), adjust = "tukey")

##===============================================
## How much of the variance is explained by the fixed and random factors? ----
# Requires library(ape)

(v1 <- varcomp(Model.3, TRUE, FALSE)*100)  ## in percent (sums up to 100%)
(v2 <- varcomp(Model.3, FALSE, FALSE))

##===============================================
## Data summary ----

# Average somatic growth rate for each clone in each treatment

avg.clone <- myfile %>% 
        group_by(Genotype,PUFA) %>% 
        summarize(Mean.clone = mean(Growth_rate, na.rm = T),
                  SD.clone = sd(Growth_rate, na.rm = T))

avg.clone

# Average somatic growth rate of the entire population in each treatment

avg.pop <- avg.clone %>% 
        group_by(PUFA) %>% 
        summarize(Mean.pop = mean(Mean.clone, na.rm = T),
                  SD.pop = sd(Mean.clone, na.rm = T))

avg.pop