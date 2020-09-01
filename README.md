# Intrapopulation-variability-in-a-functional-trait

# General Information

This repository contains raw data and R code used within the following publication:

**Intrapopulation variability in a functional trait: susceptibility of Daphnia to limitation by dietary fatty acids**

Maja Ilić1,2*, Mathilde Cordellier3† and Patrick Fink1,4,5°

1 University of Cologne, Institute for Zoology, Zülpicher Strasse 47b, 50674 Köln, Germany 
2 Present address: Queen’s University Belfast, School of Biological Sciences, 19 Chlorine Gardens, Belfast BT9 5DL, United Kingdom 
3 Universität Hamburg, Institute of Zoology, Martin-Luther-King Platz 3, 20146 Hamburg, Germany
4 Helmholtz Centre of Environmental Research – UFZ, Department River Ecology, Brückstrasse 3a, 39114 Magdeburg, Germany 
5 Helmholtz Centre of Environmental Research – UFZ, Department Aquatic Ecosystem Analysis and Management, Brückstrasse 3a, 39114 Magdeburg, Germany

* Author for correspondence: 
E-mail: *maja.ilic.bio@gmail.com*
ORCID iD: 0000-0002-8387-9932

† E-mail: mathilde.cordellier@uni-hamburg.de
ORCID iD: 0000-0001-7376-4560

° E-mail: patrick.fink@ufz.de
ORCID iD: 0000-0002-5927-8977

# Raw data

Raw data can be found in the Excel file *Raw data Ilic, Cordellier and Fink, Freshwater Biology.xlsx*
Sheet *SGR*: somatic growth rate for each genotype 
Sheet *Susceptibility*: susceptibility for each genotype

**Explanation of all variables (columns):**

*Day*: Gives a day of the year 2017 on which each individual experiment was started within the experimental period (March - November 2017)                        
       Day 1: 1st of January, 2017                                                   
                                                                                    
*Experiment_Nr*: Number of the experiment for every of the genotypes (at least 3 experiments per genotype)                               
                                                                                    
*Genotype*: Laboratory code/name for each *D. longispina genotype* isolated from the Lake Klostersee (KL) and used in this study                          
                                                                                    
*Species*: Species, identified following Rusek et al. (2015).                        
           Note that 11 out of 12 genotypes are *D. longispina*, while only the genotype KL8 is a *D. longispina × cucullata* hybrid                    
                                                                                    
*PUFA*: Liposomes used to supplement the diet (green alga Acutodesmus obliquus) with following PUFAs: 
        C - control liposomes, PUFA-free;                      
        ALA - alpha-linolenic acid;                            
        EPA - eicosapentaenoic acid;                           
        ARA - arachidonic acid                                 
                             
*Growth_rate*: Somatic growth rate of *Daphnia* individuals (after 6 days of experiment) in every jar per treatment (see Liposome), calculated with the Eq. 1 (see M and M)

*Susceptibility*: Susceptibility, calculated with the Eq. 2, given in %                                                       
                  Negative value indicates that the somatic growth rate of *Daphnia* was higher in presence of the corresponding PUFA compared to the control treatment.
                  
# R code

R code is subdivided into four scripts, labelled with 01-04.

Part 1: Somatic growth rate - Linear mixed-effects models (for entire population)  
Part 2: Somatic growth rate - Linear mixed-effects models (for each clone individually)  
Part 3: Susceptibility - Linear mixed-effects models (for entire population)  
Part 4: Susceptibility - Linear mixed-effects models (for each clone individually)  

All required packages are given at the very beginning of each script.

# Final notes

For all remarks, suggestions or questions, please contact the corresponding author (me) under maja.ilic.bio@gmail.com



