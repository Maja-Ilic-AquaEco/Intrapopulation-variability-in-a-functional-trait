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
## Part 4: Susceptibility of Daphnia as response variable (individual clones)             ##  
##                                                                                        ##
## Containts code for data import and statistical analyses                                ##
##                                                                                        ##    
############################################################################################

##===============================================
## Packages ----

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(nlme)
library(broom.mixed)
library(purrr)
library(car)
library(multcomp)

##===============================================
## Set working directory ----

setwd("~/My Documents/Paper Clonal variability in susceptibility/Data and scripts following publication")

##===============================================
## Import raw data ----

df <- read_excel("Raw data Ilic, Cordellier and Fink, Freshwater Biology.xlsx", sheet = "Susceptibility")
head(df)
str(df)

##===============================================
## LME for each genotype using tidyr, broom.mixed and purrr ----

# Let control treatment be the baseline

df$PUFA <- factor(df$PUFA, levels = c("C","ALA","EPA","ARA"))

# Fit a linear mixed-effects model for each genotype with Susceptibility as response variable, 
# PUFA as fixed effect and Day of the year as random effect

lme_Susceptibility <- df %>%
  dplyr::select(Day,Genotype,PUFA,Susceptibility) %>% 
  nest(-Genotype) %>% 
  mutate(
    fit = map(data, ~ lme(Susceptibility ~ PUFA,
                          random = ~ 1|Day,
                          data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  ) %>% 
  unnest(tidied)

# Perform anova on each fit

anova_Susceptibility <- df %>%
  dplyr::select(Day,Genotype,PUFA,Susceptibility) %>% 
  nest(-Genotype) %>% 
  mutate(
    fit = map(data, ~ lme(Susceptibility ~ PUFA,
                          random = ~ 1|Day,
                          data = .x)),
    anova = map(fit, anova)
  ) %>% 
  unnest(anova)

# Note that each second line refers to the fixed factor PUFA (numDF = 3)

anova_Susceptibility %>%  filter(numDF == 3)

##===============================================
## Source functions for histogramm and QQ plot ----

source("Function ggHistNorm.R")
source("Function ggQQplot.R")

##===============================================
## Model asumptions, validation and Dunnett's posthoc test ----

dir.create("ggHistNorm Residuals Susceptibility")
dir.create("ggQQplot Residuals Susceptibility")
dir.create("Variance homogeneity Residuals Susceptibility")
dir.create("Homoscedasticity for Susceptibility")

genotypes <- unique(df$Genotype)

anova.F <- c()
anova.p <- c()
shapiro.p <- c()
levene.p <- c()

for(i in 1:length(genotypes)){
  
  mygenotype <- df %>% 
    dplyr::select(Day,Genotype,PUFA,Susceptibility) %>% 
    filter(Genotype == genotypes[i]) %>% 
    droplevels()
  
  mygenotype$PUFA <- factor(mygenotype$PUFA, levels = c("C","ALA","EPA","ARA"))
  
  lme1 <- lme(Susceptibility ~ PUFA,
              random = ~ 1|Day,
              data = mygenotype)
  
  # Check normality of the residuals
  
  mygenotype$Residuals <- residuals(lme1, type = "pearson")
  mygenotype$Fitted <- fitted(lme1)
  
  shapiro.p[i] <- shapiro.test(mygenotype$Residuals)$p.value
  
  ggHistNorm.resid <- ggHistNorm(mygenotype,mygenotype$Residuals,paste("Standardized Residuals for Susceptibility - Genotype",genotypes[i]),0.5)[[1]]
  
  ggsave(paste0("ggHistNorm Residuals Susceptibility/Genotype ",genotypes[i],".png"),
         ggHistNorm.resid, width = 10, height = 8)
  
  ggQQplot.resid <- ggQQplot(mygenotype,mygenotype$Residuals,paste("Standardized Residuals for Susceptibility - Genotype",genotypes[i]))
  
  ggsave(paste0("ggQQplot Residuals Susceptibility/Genotype ",genotypes[i],".png"),
         ggQQplot.resid, width = 6, height = 6)
  
  # Check homogeneity of variances of residuals
  
  levene.p[i] <- leveneTest(Residuals ~ PUFA, data = mygenotype)$"Pr(>F)"[1]
  
  boxplot.liposome <- ggplot(mygenotype, aes(x = PUFA, y = Residuals, fill = PUFA, color = PUFA)) +
    geom_boxplot(alpha = 0.5) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(y = "Standardized Residuals",
         title = paste("Genotype",genotypes[i]),
         subtitle = paste("Levene's test p-value:",round(levene.p[i], digits = 3)))
  
  ggsave(paste0("Variance homogeneity Residuals Susceptibility/Genotype ",genotypes[i],".png"),
         boxplot.liposome, width = 6, height = 6)
  
  # Check homoscedasticity
  
  homoscedasticity.plot <- ggplot(mygenotype, aes(x = Fitted, y = Residuals, fill = PUFA, color = PUFA)) +
    geom_point(size = 3, alpha = 0.5) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Fitted values",
         y = "Standardized Residuals",
         title = paste("Genotype",genotypes[i]))
  
  ggsave(paste0("Homoscedasticity for Susceptibility/Genotype ",genotypes[i],".png"),
         homoscedasticity.plot, width = 6, height = 6)
  
  # Run anova
  
  aov1 <- anova(lme1)
  
  anova.F[i] <- aov1$`F-value`[2]
  anova.p[i] <- aov1$`p-value`[2]
  
  # Perform Dunnett's test on each fit
  
  dunnett <- summary(glht(lme1, linfct = mcp(PUFA = "Dunnett")))
  
  results.posthoc <- data.frame("Genotype" = rep(genotypes[i],3),
                                "Comparison" = row.names(data.frame(dunnett$linfct)),
                                "Estimate" = dunnett$test$coefficients,
                                "Std.Error" = dunnett$test$sigma,
                                "t.value" = dunnett$test$tstat,
                                "p.value" = dunnett$test$pvalue)
  
  if(i == 1){
    results.final <- results.posthoc
  }
  
  if(i > 1){
    results.final <- rbind(results.final,results.posthoc)
  }
  
  rm(mygenotype,lme1,aov1,dunnett,results.posthoc)
}

##===============================================
## Export results ----

results.df <- cbind(as.character(genotypes),anova.F,anova.p,shapiro.p,levene.p)

write.table(results.df, "LME model asumptions and output for each genotype - Susceptibility.csv", sep = ",", row.names = F)

write.table(results.final, "Dunnett posthoc test for each genotype - Susceptibility.csv", sep = ",", row.names = F)