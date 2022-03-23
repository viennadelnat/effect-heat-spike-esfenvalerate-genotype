#How a heat spike shapes pesticide toxicity differs among clones of the water flea Daphnia magna
#Vienna Delnat, Julie Verheyen, Ine Van Hileghem and Robby Stoks
#Journal (Year)
#R code tested on 23/03/2022


#####Packages#####

install.packages("car")
install.packages("lme4")
install.packages("afex")
install.packages("emmeans")
install.packages("effects")
install.packages("ggpubr")
install.packages("tidyverse")
install.packages("rlang")
install.packages("pbkrtest")

library(car)     
library(lme4)
library(afex)
library(emmeans)
library(effects)     
library(ggpubr)
library(tidyverse)
library(pbkrtest)

sessionInfo()


#####Datasets#####

#Intrinisc population growth rate calculations using the Euler Lotka equation 
#Note that the Number of Days per ID (jar) differ as the jars are kept for the number of days till the mother released second brood
dataPopGrowth=read.csv("./Delnat-et-al_EulerLotka.csv", sep=",", na.strings=c(""))
dataPopGrowth$Pesticide=factor(dataPopGrowth$Pesticide, levels=c("solvent","esfen"))
Factors <- c("Clone", "HeatSpike", "StartDate")
dataPopGrowth[Factors] <- do.call(cbind.data.frame, lapply(dataPopGrowth[Factors], as.factor))
dataPopGrowth$Day = as.numeric(as.character(dataPopGrowth$Day))
str(dataPopGrowth)

#AllData (averages per jar)
data=read.csv("./Delnat-et-al_MeanData_PopulationGrowthRate_TimeToMaturation_GrazingRate.csv", sep=",", na.strings=c(""))
#set variables as factors
data$Pesticide=factor(data$Pesticide, levels=c("solvent","esfen"))
Factors <- c("Clone", "HeatSpike", "StartDate")
data[Factors] <- do.call(cbind.data.frame, lapply(data[Factors], as.factor))
str(data)

#Survival individual measurements (binomial data)
dataSurv=read.csv("./Delnat-et-al_IndividualData_Survival.csv", sep=",", na.strings=c(""))
dataSurv$Pesticide=factor(dataSurv$Pesticide, levels=c("solvent","esfen"))
Factors <- c("Clone", "HeatSpike", "StartDate", "ID")
dataSurv[Factors] <- do.call(cbind.data.frame, lapply(dataSurv[Factors], as.factor))
str(dataSurv)
dataSurvMat=subset(dataSurv, SurvivalMaturationFirst!="")
str(dataSurvMat)

#Body size individual measurements
dataBodySize=read.csv("./Delnat-et-al_IndividualData_BodySize.csv", sep=",", na.strings=c(""))
dataBodySize$Pesticide=factor(dataBodySize$Pesticide, levels=c("solvent","esfen"))
Factors <- c("Clone", "HeatSpike", "StartDate", "ID")
dataBodySize[Factors] <- do.call(cbind.data.frame, lapply(dataBodySize[Factors], as.factor))
str(dataBodySize)

#CTmax individual measurements
dataCTmax=read.csv("./Delnat-et-al_IndividualData_CTmax.csv", sep=",", na.strings=c(""))
dataCTmax$Pesticide = factor(dataCTmax$Pesticide, levels=c("solvent","esfen"))
Factors <- c("Clone", "HeatSpike", "StartDate", "ID", "NoBodySize")
dataCTmax[Factors] <- do.call(cbind.data.frame, lapply(dataCTmax[Factors], as.factor))
#Only use CTmax values for individuals for which the body size is available for correction (ca. 6% data lost)
dataCTmax=subset(dataCTmax,NoBodySize=="no")
str(dataCTmax)

#MDR data
dataMDR=read.csv("./Delnat-et-al_MDR.csv", sep=",", na.strings=c(""))
dataMDR["Clone"]=do.call(cbind.data.frame, lapply(dataMDR["Clone"], as.factor))
str(dataMDR)


######Euler Lotka calculations for intrinsic population growth rate (r)######

#remark - jars (IDs) with on a certain day no mothers alive will cause an error (is an NA, which gives an error with the uniroot), 
#if your mothers died before the end of the second brood, discard these jars 
#or use these jars till the last day that they you had living mothers (can still be discarded later if needed)

# calculate lx and mx are needed to calculate rm
daphnia_growth <- dataPopGrowth %>%
  group_by(Clone, HeatSpike, Pesticide, ID) %>%
  mutate(lx = NumberFemales/max(NumberFemales)) %>%
  mutate(mx = (NumberBrood/NumberFemales)) %>%
  ungroup()
daphnia_growth=as.data.frame(daphnia_growth)
str(daphnia_growth)

#new data frame of all jars present in data frame
Jar=unique(daphnia_growth["ID"])
#list of unique jar IDs
jarList=Jar$ID
#number of unique jars
NumberOfJars=length(jarList)
#list of column names
columns=c("ID", "ReproductivePerformance", "Lambda")
#number of unique columns
NumberOfColumns=length(columns)
#empty matrix with correct number of columns and rows
RPdata=matrix(ncol = NumberOfColumns, nrow = NumberOfJars)
#add column names to new empty matrix
colnames(RPdata)=columns
#for loop to calculate the reproductive performance using the Euler-Lotka equation
for (i in seq(1,NumberOfJars,1)){
  SingleJar = jarList[i]
  #subset per jar
  daphnia_growthFor = subset(daphnia_growth, ID == SingleJar)
  # Calculate rm (intrinsic growth rate) using Euler-Lotka equation 
  x <- c(daphnia_growthFor$Day) 
  L <- c(daphnia_growthFor$lx) 
  m <- c(daphnia_growthFor$mx) 
  # a vector containing the end-points of the interval to be searched for the root
  r.range<- c(0, 5) 
  # introduce the Euler-Lotka equation for solution (root) 
  eulerlotka <- function(r) sum(L * m * exp(-r*x)) - 1 
  # define uniroot and set acuracy of estimation (tol)
  res <- uniroot(f = eulerlotka, interval = r.range, extendInt = "yes", tol = 1e-8) 
  #add jar ID to new matrix
  RPdata[i, "ID"]=SingleJar
  #extract and round the reproductive performance
  RPdata[i, "ReproductivePerformance"]=round(res$root, digits = 2)
  #calculate Lambda (exp(RP))
  RPdata[i, "Lambda"]=round(exp(res$root), digits = 2)
}

#Set reproductive performance as data frame and set correct data types for each variable
RPdataDF=as.data.frame(RPdata)

#Left join to add Clone, HeatSpike, Pesticide and StartDate to the RPdataDF
daphnia_growth_subsetcolumns=daphnia_growth[,c(1,5:8)]
daphnia_growth_subsetcolumnsUnique=distinct(daphnia_growth_subsetcolumns)
RPdataDF=left_join(RPdataDF, daphnia_growth_subsetcolumnsUnique, by = "ID")
str(RPdataDF)

#Manually enter RP = 0 for jars were no females survived, so do not use this dataframe to do your statistics on unless all jars had at least a single surviving female
#Convert StartDate to character, as one of the dates is a new level
RPdataDF$StartDate=as.character(RPdataDF$StartDate)
#Add new rows to the dataframe
RPdataDF[nrow(RPdataDF) + 1,] <- c(48,"",0,"M7","20", "esfen","10/10/2020")
RPdataDF[nrow(RPdataDF) + 1,] <- c(139,"",0,"M8","20", "esfen","19/10/2020")
#Convert StartDate back to factor
RPdataDF$StartDate=as.factor(RPdataDF$StartDate)

#Set correct data type
Numerics <- c("ReproductivePerformance", "Lambda")
RPdataDF[Numerics] <- do.call(cbind.data.frame, lapply(RPdataDF[Numerics], as.numeric))
str(RPdataDF)

#Save dataframe as csv file in working directory on pc
write.csv(RPdataDF, "Delnat-et-al_ReproductivePerformance.csv")


######Random factors and covariates######

#StartDate was added as random factor for all endpoints as the vials were started across multiple days

#ID (=vial) was added to the models as a random factor to take into account that Daphnia from the same vial were not independent 
#as individuals were used as the unit of replication for the survival endpoints, CTmax and body size at maturation.

#CTmax and grazing rate were corrected for the body size after the measurement of CTmax (5-day old juvenile Daphnia) 
#or of grazing rate (24-48h mature Daphnia) respectively by adding it as a covariate. 


######Survival 48 h after the start of the pesticide pulse######

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
glmerSurv48=glmer(SurvivalPesticide48h ~ HeatSpike*Pesticide*Clone + (1|StartDate) + (1|ID), data=dataSurv, na.action=na.omit, family=binomial(link=logit),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerSurv48, type=3)
#Assumption - Dispersion parameter
glmSurv48=glm(SurvivalPesticide48h ~ HeatSpike*Pesticide*Clone, data=dataSurv, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmSurv48) 
#Effects single stressor: Pesticide - Pesticide:Clone not significant - Contrasts based on 2-way interaction (HeatSpike:Pesticide) 
interactS1<-pairs(emmeans(glmerSurv48, ~HeatSpike|Pesticide, adjust="none"))
interactS2<-pairs(emmeans(glmerSurv48, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactS1,interactS2, adjust="fdr")
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone)
contrastS<-pairs(emmeans(glmerSurv48, ~HeatSpike*Pesticide|Clone, adjust="none"))
rbind(contrastS, adjust="fdr")


######Reproductive performance (Lambda=exp(RP))######

#Lambda = exp(RP) (used when some vials have no survivors till first brood)

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
lmLambda=lmer(Lambda ~ HeatSpike*Pesticide*Clone + (1|StartDate), data=RPdataDF, na.action=na.omit)
#Assumption of normality is not ok (W=0.72) --> Box-Cox transformation needed
boxcox=lm(Lambda ~HeatSpike*Pesticide*Clone, data=RPdataDF, na.action=na.omit) 
#Plot profile Log-likelihood: range determined by 'seq' --> change the range till you find the peak
car::boxCox(boxcox, family="yjPower", plotit=TRUE, lambda = seq(4.5,5.3 , length = 10))
#Plot profile Log-likelihood: peak = 4.9 <-- power used in box cox transformation
bcLambda <- car::yjPower(RPdataDF$Lambda, 4.9)
#Add Box-Cox transformed Lambda to your dataset
RPdataDF$bcLambda=bcLambda
#Run model with Box-Cox transformed Lambda
lmLambda=lmer(bcLambda ~ HeatSpike*Pesticide*Clone + (1|StartDate), data=RPdataDF, na.action=na.omit)
Anova(lmLambda, type=3)
#Assumptions
shapiro.test(resid(lmLambda))                    
hist(resid(lmLambda))      
leveneTest(Lambda ~ HeatSpike*Pesticide*Clone, data = RPdataDF) 
leveneTest(bcLambda ~ HeatSpike*Pesticide*Clone, data = RPdataDF) 
#influential observations
cd=cooks.distance(lmLambda)
inflobs=which(cd>1);inflobs 
influenceIndexPlot(lmLambda,vars="Cook") 
#Emmeans for Figure 2A in main text using Excel and PowerPoint
emmeans(lmLambda, ~Clone*HeatSpike*Pesticide, adjust="fdr")
#Effects single stressor: Pesticide - Pesticide:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone) 
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone)
contrastL<-pairs(emmeans(lmLambda, ~HeatSpike*Pesticide|Clone, adjust="none"))
rbind(contrastL, adjust="fdr")
#Interaction effect between both stressors - Clone independent
interactL1<-pairs(emmeans(lmLambda, ~Pesticide|HeatSpike, adjust="none"))
interactL2<-pairs(emmeans(lmLambda, ~HeatSpike|Pesticide, adjust="none"))
rbind(interactL1,interactL2, adjust="fdr")

#Plot (not used in manuscript)
LambdaData <- summary(emmeans(lmLambda, ~ HeatSpike*Pesticide*Clone, type = "response"))
png("Figure_Lambda.png", width = 18, height = 18, unit = "cm", res=300)
LambdaPlot <- ggplot(LambdaData, aes(x = Clone, y = emmean, col = Pesticide, shape = HeatSpike)) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5)) +
  #Set the position and size of the points manually
  geom_point(size = 5, position = position_dodge(.5)) +
  #set shape and label names
  scale_shape(name = "", labels = c("No heat spike", "Heat spike")) +
  #set colors and label names
  scale_color_manual(values=c("#0044FF", "#FFB700"), name = "", labels = c("Solvent control", "Esfenvalerate")) +
  #Set axis tick labels
  scale_x_discrete(breaks=c("B7", "M4", "M5","M6", "M7", "M8"),
                   labels=c("1", "2", "3", "4", "5", "6")) +
  #Set y-label
  ylab(expression("Reproductive performance (box cox, lambda=exp(RP))")) +
  #Set x-label
  xlab(expression("Daphnia clone (genotype)")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "bottom",
    #Set font size, color and margins of the legend title
    legend.title=element_text(size=20, colour = "grey29"),
    #Set font size, color and margins of the legend
    legend.text = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "grey29"),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "grey29"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "grey29"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt"), colour = "grey29"),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt"), colour = "grey29"),
    #Set colour of ticks
    axis.ticks = element_line(colour = "grey29"),
    #Set colour of axis line
    axis.line = element_line(colour = "grey29"),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
LambdaPlot
dev.off()


######Time to maturation######

#Time to maturation indicates the time till half of the surviving Daphnia in a vial were mature

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
lmTimeMatHalf=lmer(TimeMaturationHalfSurvival~ HeatSpike*Pesticide*Clone + (1|StartDate), data=data, na.action=na.omit) 
Anova(lmTimeMatHalf, type=3)
#Assumptions
shapiro.test(resid(lmTimeMatHalf))                    
hist(resid(lmTimeMatHalf))      
leveneTest(TimeMaturationHalfSurvival ~ HeatSpike*Pesticide*Clone, data = data) 
#influential observations
cd=cooks.distance(lmTimeMatHalf)
inflobs=which(cd>1);inflobs 
influenceIndexPlot(lmTimeMatHalf,vars="Cook") 
#Emmeans for Figure 2B in main text using Excel and PowerPoint
emmeans(lmTimeMatHalf, ~Clone*HeatSpike*Pesticide, adjust="fdr")
#Effects single stressor: Pesticide - Pesticide:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone) 
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone)
contrastT<-pairs(emmeans(lmTimeMatHalf, ~HeatSpike*Pesticide|Clone, adjust="none"))
rbind(contrastT, adjust="fdr")

#Determination of the interaction type per clone with separate GLMERs per clone as HS:Pesticide:Clone is significant in the main model
#Analyzing the interaction type between a heat spike and esfenvalerate exposure of Clone 1
data1=subset(data, Clone=="B7") #Clone 1
lmTimeMatHalf=lmer(TimeMaturationHalfSurvival~ HeatSpike*Pesticide + (1|StartDate), data=data1, na.action=na.omit) 
Anova(lmTimeMatHalf, type=3)
interactT1a<-pairs(emmeans(lmTimeMatHalf, ~HeatSpike|Pesticide, adjust="none"))
interactT1b<-pairs(emmeans(lmTimeMatHalf, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactT1a,interactT1b, adjust="fdr")
#Analyzing the interaction type between a heat spike and esfenvalerate exposure of Clone 2
data2=subset(data, Clone=="M4") #Clone 2
lmTimeMatHalf=lmer(TimeMaturationHalfSurvival~ HeatSpike*Pesticide + (1|StartDate), data=data2, na.action=na.omit) 
Anova(lmTimeMatHalf, type=3)
#Analyzing the interaction type between a heat spike and esfenvalerate exposure of Clone 3
data3=subset(data, Clone=="M5") #Clone 3
lmTimeMatHalf=lmer(TimeMaturationHalfSurvival~ HeatSpike*Pesticide + (1|StartDate), data=data3, na.action=na.omit) 
Anova(lmTimeMatHalf, type=3)
#Analyzing the interaction type between a heat spike and esfenvalerate exposure of Clone 4
data4=subset(data, Clone=="M6") #Clone 4
lmTimeMatHalf=lmer(TimeMaturationHalfSurvival~ HeatSpike*Pesticide + (1|StartDate), data=data4, na.action=na.omit) 
Anova(lmTimeMatHalf, type=3)
plot(effect(mod=lmTimeMatHalf, term="HeatSpike*Pesticide"))
#Analyzing the interaction type between a heat spike and esfenvalerate exposure of Clone 5
data5=subset(data, Clone=="M7") #Clone 5
lmTimeMatHalf=lmer(TimeMaturationHalfSurvival~ HeatSpike*Pesticide + (1|StartDate), data=data5, na.action=na.omit) 
Anova(lmTimeMatHalf, type=3)
#Analyzing the interaction type between a heat spike and esfenvalerate exposure of Clone 6
data6=subset(data, Clone=="M8") #Clone 6
lmTimeMatHalf=lmer(TimeMaturationHalfSurvival~ HeatSpike*Pesticide + (1|StartDate), data=data6, na.action=na.omit) 
Anova(lmTimeMatHalf, type=3)
interactT6a<-pairs(emmeans(lmTimeMatHalf, ~HeatSpike|Pesticide, adjust="none"))
interactT6b<-pairs(emmeans(lmTimeMatHalf, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactT6a,interactT6b, adjust="fdr")


######Body size at maturation######

#Body size at maturation uses the body sizes measured for the grazing rate of 24-48h mature Daphnia

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
lmSize=lmer(BodySize ~ HeatSpike*Pesticide*Clone + (1|StartDate) + (1|ID), data=dataBodySize, na.action=na.omit)
Anova(lmSize, type=3, white.adjust = TRUE)
#white.adjust = TRUE since assumption of homogeneity was not met
#Assumptions
shapiro.test(resid(lmSize))                    
hist(resid(lmSize))      
leveneTest(BodySize ~ HeatSpike*Pesticide*Clone, data = dataBodySize) 
#Rule of thumb: if smallest and largest variance do not differ with more than factor 5, assumption of homogeneity is still met
aggregate(BodySize ~ HeatSpike*Pesticide*Clone, data = dataBodySize, var) 
#influential observations
cd=cooks.distance(lmSize)
inflobs=which(cd>1);inflobs 
influenceIndexPlot(lmSize,vars="Cook") 
#Emmeans for Figure 2C in main text using Excel and PowerPoint
emmeans(lmSize, ~Clone*HeatSpike*Pesticide, adjust="fdr")
#Effects single stressor: Pesticide - Pesticide:Clone not significant - Contrasts based on 2-way interaction (HeatSpike:Pesticide) 
interactB1<-pairs(emmeans(lmSize, ~HeatSpike|Pesticide, adjust="none"))
interactB2<-pairs(emmeans(lmSize, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactB1,interactB2, adjust="fdr")
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone)
contrastB<-pairs(emmeans(lmSize, ~HeatSpike*Pesticide|Clone, adjust="none"))
rbind(contrastB, adjust="fdr")


######CTmax######

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
lmCTmax=lmer(CTmax~ HeatSpike*Pesticide*Clone + BodySize + (1|StartDate) + (1|ID), data=dataCTmax, na.action=na.omit) 
Anova(lmCTmax, type=3) 
#Assumptions
shapiro.test(resid(lmCTmax))                    
hist(resid(lmCTmax))      
leveneTest(CTmax ~ HeatSpike*Pesticide, data = dataCTmax) 
#influential observations
cd=cooks.distance(lmCTmax)
inflobs=which(cd>1);inflobs 
influenceIndexPlot(lmCTmax,vars="Cook") 
#Emmeans for Figure 3A in main text using Excel and PowerPoint
emmeans(lmCTmax, ~Clone*HeatSpike*Pesticide, adjust="fdr")
#Effects single stressor: Pesticide - Pesticide:Clone not significant - Contrasts based on 2-way interaction (HeatSpike:Pesticide) 
#Effects single stressor: HeatSpike - HS:Clone not significant - Contrasts based on 2-way interaction (HeatSpike:Pesticide)
interactC1<-pairs(emmeans(lmCTmax, ~HeatSpike|Pesticide, adjust="none"))
interactC2<-pairs(emmeans(lmCTmax, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactC1,interactC2, adjust="fdr")
#Main clone effect (without Clone being in an interaction with one or both stressors)
contractC<-pairs(emmeans(lmCTmax, ~Clone, adjust="none"))
rbind(contractC, adjust="fdr")


######Grazing rate######

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
lmGrazing=lmer(GrazingRateCellsmL ~ HeatSpike*Pesticide*Clone + BodySizeGrazingRate + (1|StartDate), data=data, na.action=na.omit)
Anova(lmGrazing, type=3, white.adjust = TRUE)
#white.adjust = TRUE since assumption of homogeneity was not met
#Assumptions
shapiro.test(resid(lmGrazing))                    
hist(resid(lmGrazing))      
leveneTest(GrazingRateCellsmL ~ HeatSpike*Pesticide*Clone, data = data) 
#Rule of thumb: if smallest and largest variance do not differ with more than factor 5, assumption of homogeneity is still met
aggregate(GrazingRateCellsmL ~ HeatSpike*Pesticide*Clone, data = data, var) 
#influential observations
cd=cooks.distance(lmGrazing)
inflobs=which(cd>1);inflobs 
influenceIndexPlot(lmGrazing,vars="Cook") 
#Emmeans for Figure 3B in main text using Excel and PowerPoint
emmeans(lmGrazing, ~Clone*HeatSpike*Pesticide, adjust="fdr")
#Effects single stressor: Pesticide - Pesticide:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone) 
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone)
contrastG<-pairs(emmeans(lmGrazing, ~HeatSpike*Pesticide|Clone, adjust="none"))
rbind(contrastG, adjust="fdr")
#Effects single stressor: Pesticide - no contrast (20 solvent - 20 esfen) based on 3-way interaction (HeatSpike:Pesticide:Clone) was significant for any of the six clones
interactG1<-pairs(emmeans(lmGrazing, ~HeatSpike|Pesticide, adjust="none"))
interactG2<-pairs(emmeans(lmGrazing, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactG1,interactG2, adjust="fdr")


######Survival after the heat spike exposure (Appendix C)######

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
glmerSurvHS=glmer(SurvivalHeatSpike ~ HeatSpike*Clone + (1|StartDate) + (1|ID), data=dataSurv, na.action=na.omit, family=binomial(link=logit),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerSurvHS, type=3)
#Assumption - Dispersion parameter
glmSurvHS=glm(SurvivalHeatSpike ~ HeatSpike*Clone, data=dataSurv, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmSurvHS) 
#Emmeans for Figure C1 in Appendix C using Excel and PowerPoint
emmeans(glmerSurvHS, ~Clone*HeatSpike, adjust="fdr",type="response")
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 2-way interaction (HeatSpike:Clone)
contrastCS<-pairs(emmeans(glmerSurvHS, ~HeatSpike|Clone, adjust="none"))
rbind(contrastCS, adjust="fdr")


######Survival at maturation (Appendix C)######

#Survival at maturation indicates the survival when the first mature Daphnia in a vial was present

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
glmerSurvMat=glmer(SurvivalMaturationFirst ~ HeatSpike*Pesticide*Clone + (1|StartDate) + (1|ID), data=dataSurvMat, na.action=na.omit, family=binomial(link=logit),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerSurvMat, type=3)
#Assumption - Dispersion parameter
glmSurvMat=glm(SurvivalMaturationFirst ~ HeatSpike*Pesticide*Clone, data=dataSurvMat, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmSurvMat) 
#Emmeans for Figure C2 in Appendix C using Excel and PowerPoint
emmeans(glmerSurvMat, ~Clone*HeatSpike*Pesticide, adjust="fdr",type="response")
#Effects single stressor: Pesticide - Pesticide:Clone not significant - Contrasts based on 2-way interaction (HeatSpike:Pesticide) 
interactCM1<-pairs(emmeans(glmerSurvMat, ~HeatSpike|Pesticide, adjust="none"))
interactCM2<-pairs(emmeans(glmerSurvMat, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactCM1,interactCM2, adjust="fdr")
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone)
contrastCM<-pairs(emmeans(glmerSurvMat, ~HeatSpike*Pesticide|Clone, adjust="none"))
rbind(contrastCM, adjust="fdr")


######Fecundity per female (Appendix C)######

#Fecundity per female is the sum of the number of juveniles from the first and second brood

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()
lmFecundity=lmer(FecundityPerFemale~ HeatSpike*Pesticide*Clone + (1|StartDate), data=data, na.action=na.omit) 
#white.adjust = TRUE since assumption of homogeneity was not met
Anova(lmFecundity, type=3, white.adjust=TRUE)
#Assumptions
shapiro.test(resid(lmFecundity))                    
hist(resid(lmFecundity))      
leveneTest(FecundityPerFemale ~ HeatSpike*Pesticide*Clone, data = data) 
#Rule of thumb: if smallest and largest variance do not differ with more than factor 5, assumption of homogeneity is still met
aggregate(FecundityPerFemale ~ HeatSpike*Pesticide*Clone, data = data, var) 
#influential observations
cd=cooks.distance(lmFecundity)
inflobs=which(cd>1);inflobs 
influenceIndexPlot(lmFecundity,vars="Cook") 
#Emmeans for Figure C3 in Appendix C using Excel and PowerPoint
emmeans(lmFecundity, ~Clone*HeatSpike*Pesticide, adjust="fdr")
#Effects single stressor: Pesticide - Pesticide:Clone not significant - Contrasts based on 2-way interaction (HeatSpike:Pesticide) 
interactCF1<-pairs(emmeans(lmFecundity, ~HeatSpike|Pesticide, adjust="none"))
interactCF2<-pairs(emmeans(lmFecundity, ~Pesticide|HeatSpike, adjust="none"))
rbind(interactCF1,interactCF2, adjust="fdr")
#Effects single stressor: HeatSpike - HS:Clone significant - Contrasts based on 3-way interaction (HeatSpike:Pesticide:Clone)
contrastCF<-pairs(emmeans(lmFecundity, ~HeatSpike*Pesticide|Clone, adjust="none"))
rbind(contrastCF, adjust="fdr")


######Correlation between interaction strength and pesticide/heat tolerance######

#Interaction strength = MDR (model deviation ratio) = (predicted combined effect based on IA) / (observed combined effect)
#Pesticide tolerance = SurvivalPesticide48h = Survival 48 hours after the start of the pesticide pulse (subset: only not preceded by heat spike)
#Heat tolerance = CTmax (subset: only solvent control not preceded by heat spike)

#non-parametric correlation test
cor.test(dataMDR$MDR, dataMDR$SurvivalPesticide48h, method="spearman",exact=FALSE)
cor.test(dataMDR$MDR, dataMDR$CTmax, method="spearman",exact=FALSE)

#Figure 3
#horizontal straight line indicating additivity plus annotations 'synergism' and 'antagonism' added using PowerPoint 
ggscatter(dataMDR, x = "SurvivalPesticide48h", y = "MDR", ylim=c(0.4,1.3),
          conf.int = FALSE, cor.coef = FALSE, cor.method = "spearman", add = "reg.line",
          xlab = "Pesticide tolerance (% survival 48 hours after esfenvalerate pulse)", ylab = "Model Deviation Ratio")
ggscatter(dataMDR, x = "CTmax", y = "MDR", ylim=c(0.4,1.3),
          conf.int = FALSE, cor.coef = FALSE, cor.method = "spearman",
          xlab = "Heat tolerance (°C of CTmax under standard conditions)", ylab = "Model Deviation Ratio")


######Correlation between heat and pesticide tolerance######
cor.test(dataMDR$CTmax, dataMDR$SurvivalPesticide48h, method="spearman",exact=FALSE)

ggscatter(dataMDR, x = "CTmax", y = "SurvivalPesticide48h", 
          conf.int = FALSE, cor.coef = FALSE, cor.method = "spearman",
          xlab = "Heat tolerance (°C of CTmax under standard conditions)", ylab = "Pesticide tolerance (% survival 48 hours after esfenvalerate pulse)")

######Save Rdata######
save.image(file="Delnat et al_effect-heat-spike-esfenvalerate-genotype_20220323.Rdata")
