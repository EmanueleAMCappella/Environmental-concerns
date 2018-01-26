
# This analysis is split into two main parts:
# 1) estimation of a regression model to predict the willingness to pay for environmental concerns. I compared a model
# after multiple data imputation and a simpler listwise deletion
# 2) Factorial analysis (both EFA and CFA) on the questionnaire scale of willingness to pay

# For more details on the analyses and the data source, please read the README document



###############################################################################################
#                                 1. REGRESSION MODEL
###############################################################################################


#CHOOSING THE VARIABLES OF INTEREST: AGE, INCOME & WILLINGNESS TO PAY (WTP)
#set directory and open file 
setwd("C:/Users/EAMC/Desktop/PRINCIPLES_Coursework")
survey<- read.table("dfR.csv", header=TRUE, sep=",")
View(survey)

# Removing variables I won't use
survey<- survey[-c(1:11, 13:17, 19, 20)]
View(survey)

#eventual csv file generation
#write.csv(survey, file = "WVsurvey.csv")


# Data visualization and transformation - before regression 
#transforming dataframe columns to double/numeric object to perform the preliminary analyses
ageL<- na.omit(survey$age)
incomeL<- na.omit(survey$income)
wtpTotL<- na.omit(survey$wtpTot)

#inspect kurtosis and skewness with moments
library("moments")
skewness(ageL)
kurtosis(ageL)
skewness(wtpTotL)
kurtosis(wtpTotL)
skewness(incomeL)
kurtosis(incomeL) 

#correlation matrix
cor(survey, use = "complete.obs")

#plot age
par(mfrow=c(1,2))
h<- hist(ageL, breaks= 50, main = paste("Age Histogram & normal curve"), col= 'whitesmoke', xlab="Age")
xfit<-seq(min(ageL),max(ageL),length=40) 
yfit<-dnorm(xfit,mean=mean(ageL),sd=sd(ageL)) 
yfit <- yfit*diff(h$mids[1:2])*length(ageL) 
lines(xfit, yfit, col="slateblue1", lwd=2)

#density plot
plot(density(ageL), main="Kernel Density of Age")
polygon(density(ageL), col="whitesmoke", border="slateblue1", lwd=2)

par(mfrow=c(1,1))
#income
h1<- hist(incomeL, breaks= 10, main = paste("Income Histogram with normal curve"), col= 'gold', xlab="Income")
xfit1<-seq(min(incomeL),max(incomeL),length=50) 
yfit1<-dnorm(xfit1,mean=mean(incomeL),sd=sd(incomeL)) 
yfit1 <- yfit1*diff(h1$mids[1:2])*length(incomeL) 
lines(xfit1, yfit1, col="darkred", lwd=2)

#wtpTot
h2<- hist(wtpTotL, breaks = 12, main = paste("wtpTot Histogram & normal curve"), col= 'deepskyblue4', xlab="wtpTot")
xfit<-seq(min(wtpTotL),max(wtpTotL),length=40) 
yfit<-dnorm(xfit,mean=mean(wtpTotL),sd=sd(wtpTotL)) 
yfit <- yfit*diff(h2$mids[1:2])*length(wtpTotL) 
lines(xfit, yfit, col="gray0", lwd=2)


##############################################################################
#                   1.A. MISSING DATA ANALYSIS 
##############################################################################
survey<- read.table("dfR.csv", header=TRUE, sep=",")
View(survey)
survey<- survey[-c(1:11, 13:17, 19, 20)]
View(survey)

survey$income<- as.ordered(survey$income)
survey$wtpTot<- as.ordered(survey$wtpTot)

#Missing Data GUI for missing data visualization
MissingDataGUI(data = survey, width = 1000, height = 750)


#MULTIPLE IMPUTATION WITH PACKAGE MICE -----------------------------------------                      
# Look for missing > 5% variables
pMiss <- function(x){sum(is.na(x))/length(x)*100}

# Check each column
apply(survey,2,pMiss)

# Install mice library
library(mice)

# Missing data pattern
md.pattern(survey)

#Install VIM for graph on missing values distribution
library("VIM", lib.loc="~/R/win-library/3.4")
# Plot of missing data pattern
aggr_plot <- aggr(survey, col=c('darkslategray3','aliceblue'), numbers=TRUE, sortVars=TRUE, labels=names(survey), cex.axis=0.9, gap=2, ylab=c("Histogram of missing data","Pattern of missing data"))
#2 plots together
aggr(survey, delimiter="_imp", combined=TRUE, col=c('darkslategray3','aliceblue'), sortVars=TRUE, sortCombs=TRUE, numbers=TRUE, cex.axis=0.95, ylab=c("Histogram and Pattern of missing data combined"), prop=c(FALSE, FALSE), only.miss=TRUE)


# Impute missing data using mice
#change class of variables and the corresponding method of imputation
#method: you can define the method of imputation for each columns. If not, default take into consideration 
#the factor class that can be seen in 
str(survey)

#data imputation 
#m=5 is the default: run 5 times the imputation
#maxit = number of iteration (default = 3)

#The computation for this imputation takes a lot of time with this database...
imputeD <- mice(survey,m=5,maxit= 5, seed=333, method = c("pmm", "pmm", "pmm", "pmm"))

#details on the imputation
imputeD

#PLOTS RESULTS OF IMPUTATION -------------------------------------------------------------------------------
#plot the convergence of the algorithm
#In general, we would like the streams to intermingle and be free of any trends at the later iterations.
plot(imputeD)

# Get imputed data (for each variable)
imputeD$imp$wtpTot
imputeD$imp$income
imputeD$imp$age

# Get completed datasets (observed and imputed)
completedData1 <- complete(imputeD,1)
summary(completedData1)

# Plots  
# Density plot original vs imputed dataset
densityplot(imputeD)

#stripplot
stripplot(imputeD)

#WITH AGE NOT TRANSFORMED
# Fitting a linear model and Pooling the results 
modelfit0<- with(imputeD, lm(wtpTot~ age))
modelfit1<- with(imputeD, lm(wtpTot~ income))
modelfit2 <- with(imputeD,lm(wtpTot~ age+income))
modelfit3 <- with(imputeD,lm(wtpTot~ age+income+age*income))


#compare with simple model from dataset complete data (completedData1)
prova3<- lm(wtpTot ~ age+income+age*income, data=completedData1)
summary(prova3)
prova2<- lm(wtpTot ~ age+income, data=completedData1)
summary(prova2)

summary(modelFit1)

pool(modelfit0)
pool(modelfit1)
pool(modelfit2)
pool(modelfit3)

summary(pool(modelfit0))
summary(pool(modelfit1))
summary(pool(modelfit2))
summary(pool(modelfit3))

pool.r.squared(modelfit3, adjusted = FALSE)


##############################################################################################################
#                                     1.B. MEAN IMPUTATION    AND    LISTWISE DELETION                                                    
##############################################################################################################

#Mean imputation with R
impMR <- mice(survey, method = "mean", m = 1, maxit = 1)
fitMR0 <- with(impMR, lm(wtpTot ~ age)) 
summary(fitMR0)
fitMR1<- with(impMR, lm(wtpTot ~ income)) 
summary(fitMR1)
fitMR2<- with(impMR, lm(wtpTot ~ age+income)) 
summary(fitMR2)
fitMR3<- with(impMR, lm(wtpTot ~ age+income+age*income)) 
summary(fitMR3)


#LISTWISE DELETION
survey_listw<- read.table("dfR_listwise.csv", header=TRUE, sep=",")
View(survey_listw)

#delete useless variables
# Removing variables I won't use
survey_listw<- survey_listw[-c(1:10, 13:17, 19, 20)]  
View(survey_listw)
survey_listw<-survey_listw[-c(1)]


#--------------------------------------------------------------------------------------------------------------------
#SINCE THE REGRESSION MODEL BUILT AFTER LISTWISE DELETION IS AS ROBUST AS THE ONE COMPUTED AFTER MULTIPLE IMPUTATION, 
# I WILL USE THE LISTWISE DATA
# REGRESSION MODEL: WTP= AGE*INCOME

#fitting models regression model on listwise data
fit_list3<- lm(wtpTot ~ age+income+age*income, data=survey_listw)
summary(fit_list3)
fit_list2<- lm(wtpTot ~ age+income, data=survey_listw)
summary(fit_list2)
fit_list1<- lm(wtpTot ~ income, data=survey_listw)
summary(fit_list1)
fit_list0<- lm(wtpTot ~ age, data=survey_listw)
summary(fit_list0)

anova(fit_list0, fit_list3)
anova(fit_list2, fit_list3)


#GRAPH MODERATION ---------------------------------------------------------------------------------------------------
#split based on mean and sd of the income distribution (AGE NOT TRANSFORMED, listwise database)
x <- survey_listw$income
survey_listw$income_3g <-  case_when(x > mean(x)+sd(x) ~ "high",
                                     x < mean(x)+sd(x) & x > mean(x)-sd(x) ~ "average",
                                     x < mean(x)-sd(x) ~ "low")
count(survey_listw,  income_3g)

#graph ggplot2. Moderator: income categories. Age = x continous var, y=wtp
survey_listw %>% 
ggplot() +
aes(x = age, y = wtpTot, group = income_3g, color = income_3g) +
ggtitle("Age*income interaction") +
labs(y= "Willingness to Pay", x = "Age") +
labs(colour = "Income Categories") +
geom_point(color = "grey", alpha = .7) +
geom_smooth(method = "lm")



# INSPECT ASSUMPTION OF THE REGRESSION MODEL ------------------------------------------------------------------------------------

#for age 
install.packages("car")
library('car')
# Normality distribution of residuals, that is independence of errors. The model is:
#fit_list3<- lm(wtpTot ~ age+income+age*income, data=survey_listw)
par(mfrow=c(2,1))
hist(resid(fit_list3), breaks = 30, main='Residuals - regression model', col= 'darkseagreen2', xlab="Residuals", ylab='Frequency')
#qqplot
qqnorm(resid(fit_list3))
qqline(quantile(rnorm(1000), col= 'green')) 
durbinWatsonTest(fit_list3)

#multicollinearity
cor.test(survey_listw$age, survey_listw$income, survey_listw$wtpTot, method= 'pearson') 
vif(fit_list3) 

#gvlma package
gvmodel <- gvlma(fit_list3) 
summary(gvmodel)


##############################################################################################################
#                                         FACTORIAL ANALYSES
##############################################################################################################

#CREATE NEW SURVEY, with only the variables of interest
survey<- read.table("dfR.csv", header=TRUE, sep=",")
survey2<-na.omit(survey[-c(1, 11:24)])
#Reverse items (I forgot to reverse item wtp3)
survey2$wtp3 <- 4-survey2$wtp3
View(survey2)

#delete wtp3
#survey2<- survey2[-c(9)]

#1. Investigate discriminant and convergent validity of willingness to pay, and environmental concern 
cor(survey2)


#EFA -----------------------------------------------------------------------------------------------------
#a Estimate a plausible one, two and three factor model on the entire 
#dataset and on the three countries selected 

f1<-factanal(survey2, factors=3, rotation="promax", scores="regression")
f1
f2<-factanal(survey2, factors=3, rotation="none", scores="regression")
f2
f3<-factanal(survey2, factors=3, rotation="varimax", scores="regression")
f3

#scree plot & parallel analysis
library(nFactors)
ev <- eigen(cor(survey2)) # get eigenvalues    
ap <- parallel(subject=nrow(survey2),var=ncol(survey2), rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

#cronbach alpha for reliability check
psych::alpha(survey2, check.keys=TRUE)

#Use lavaan for Confirmatory Factor Analysis CFA -----------------------------------------------------------------
#install.packages("lavaan")
library("lavaan")

#CFA 1 factor
mod1 <- 'F1=~NA*envwrld1 + envwrld2 + envwrld3 + envcom1 + envcom2 + envcom3+ wtp1+ wtp2+ wtp3 
F1~~1*F1'
fit1<-cfa(mod1, data=survey2) 
summary(fit1, fit.measures=T)

#CFA 2 factors
model2 <- 'F1=~NA*wtp1+ wtp2+ wtp3 
F2 =~ NA*envwrld1 + envwrld2 + envwrld3 + envcom1 + envcom2 + envcom3 
F1~~1*F1
F2~~1*F2'
fit2<-cfa(model2, data=survey2)
summary(fit2, fit.measures=T)

#CFA 3 factors
model3 <- 'F1=~NA*wtp1+ wtp2 + wtp3
F2 =~ NA*envwrld1 + envwrld2 + envwrld3  
F3 =~ NA*envcom1 + envcom2 + envcom3 
F1~~1*F1
F2~~1*F2
F3~~1*F3'
fit3<-cfa(model3, data=survey2)
summary(fit3, fit.measures=T)

#CFA 4 factors
model4<- 'F1=~NA*wtp1+ wtp2
F2 =~ NA*envwrld1 + envwrld2 + envwrld3  
F3 =~ NA*envcom1 + envcom2 + envcom3 
F4 =~ NA*wtp3
F1~~1*F1
F2~~1*F2
F3~~1*F3
F4~~1*F4'
fit4<-cfa(model4, data=survey2)
summary(fit4, fit.measures=T)



###############################################################################################################
#  Compare the model fit of the best factor model for different continents & different countries
###############################################################################################################


#CREATE COUNTRY CODES WITH countrycode R Package 

survey<- read.table("dfR.csv", header=TRUE, sep=",")
#file choose: LungCapData <- read.table(file.choose(), header=T, sep="\t") 
View(survey)

#try to convert the country codes to actual names, with country code package
#install.packages('countrycode')
library("countrycode", lib.loc="~/R/win-library/3.4")

#Usage
#countrycode(sourcevar, origin, destination, warn = TRUE, custom_dict = NULL,
#            custom_match = NULL, origin_regex = FALSE)

nameCountry<- survey$country
class(nameCountry)

countryName<- countrycode(nameCountry, "wvs", "country.name.en", warn = TRUE, custom_dict = NULL,
                          custom_match = NULL, origin_regex = FALSE)

continentName<- countrycode(nameCountry, "wvs", "continent", warn = TRUE, custom_dict = NULL,
                            custom_match = NULL, origin_regex = FALSE)

#append countryName and continentName to survey dataframe
survey["countryName"]<-countryName
survey["continentName"]<- continentName

#create a new column with numbers for each continent
survey$continentNum <- NA
# Copy the data from the existing column into the new one.
survey$continentNum <- survey$continentName
#now change the values of the new column  
survey$continentNum[survey$continentNum=='Europe'] <- 1
survey$continentNum[survey$continentNum=='Africa'] <- 2
survey$continentNum[survey$continentNum=='Americas'] <- 3
survey$continentNum[survey$continentNum=='Oceania'] <- 4
survey$continentNum[survey$continentNum=='Asia'] <- 5
View(survey)

#omit na
survey<- na.omit(survey)
survey3<-survey[-c(1, 11:22)] 
View(survey3)


#c Test the best fitting model separately for Asia vs Europe ---------------------------------------------------------------------------

#create Europe subset
Europe<- survey3[survey3$continentName=='Europe', ]
#delete Europe column
Europe<- Europe[-c(10, 11)] 
View(Europe)

modelE <- 'F1=~NA*wtp1+ wtp2+ wtp3
F2 =~ NA*envwrld1 + envwrld2 + envwrld3  
F3 =~ NA*envcom1 + envcom2 + envcom3 
F1~~1*F1
F2~~1*F2
F3~~1*F3'
fita<-cfa(modelE, data=Europe)
summary(fita, fit.measures=T)

#create Asia subset
Asia<- survey3[survey3$continentName=='Asia', ]
#delete Europe column
Asia<- Asia[-c(10, 11)] 
View(Asia)

modelA <- 'F1=~NA*wtp1+ wtp2+ wtp3
F2 =~ NA*envwrld1 + envwrld2 + envwrld3         
F3 =~ NA*envcom1 + envcom2 + envcom3  
F1~~1*F1
F2~~1*F2
F3~~1*F3'
fitb<-cfa(modelA, data=Asia)
summary(fitb, fit.measures=T)

#test the fit in independent countries ------------------------------------------------------------------------------------------------

#country: norway
Norway<- survey[survey$countryName=='Norway', ]
Norway<- Norway[-c(1, 11:24)]
View(Norway)

modelN <- 'F1=~NA*wtp1+ wtp2+ wtp3
F2 =~ NA*envwrld1 + envwrld2 + envwrld3         
F3 =~ NA*envcom1 + envcom2 + envcom3  
F1~~1*F1
F2~~1*F2
F3~~1*F3'
fitN<-cfa(modelN, data=Norway)
summary(fitN, fit.measures=T)

#country: ethiopia
Ethiopia<- survey[survey$countryName=='Ethiopia', ]
Ethiopia<- Ethiopia[-c(1, 11:24)]

modelET <- 'F1=~NA*wtp1+ wtp2+ wtp3
F2 =~ NA*envwrld1 + envwrld2 + envwrld3         
F3 =~ NA*envcom1 + envcom2 + envcom3  
F1~~1*F1
F2~~1*F2
F3~~1*F3'
fitET<-cfa(modelET, data=Ethiopia)
summary(fitET, fit.measures=T)


#country: iNDIA
India<-  survey[survey$countryName=='India', ]
India<- India[-c(1, 11:24)]

modelC <- 'F1=~NA*wtp1+ wtp2+ wtp3
F2 =~ NA*envwrld1 + envwrld2 + envwrld3         
F3 =~ NA*envcom1 + envcom2 + envcom3  
F1~~1*F1
F2~~1*F2
F3~~1*F3'
fitC<-cfa(modelC, data=India)
summary(fitC, fit.measures=T)

#Country: Mexico
Mexico<-  survey[survey$countryName=='Mexico', ]
Mexico<- Mexico[-c(1, 11:24)]

modelME <- 'F1=~NA*wtp1+ wtp2+ wtp3
F2 =~ NA*envwrld1 + envwrld2 + envwrld3         
F3 =~ NA*envcom1 + envcom2 + envcom3  
F1~~1*F1
F2~~1*F2
F3~~1*F3'
fitME<-cfa(modelME, data=Mexico)
summary(fitME, fit.measures=T)



#
#######################################################################################################################
#                               OTHER ANALYSES: ANOVA on wtp FOR CONTINENTS
######################################################################################################################


#with ggplot: jitter graph
survey$jit <- jitter(survey$wtpTot, factor=3)
survey$contjit <- jitter(survey$continentNum, factor=2)
ggplot(survey, aes(y=jit, x=contjit, color=gender))+ geom_point()

#trasform continentNum e Gender in factors
survey$continentNum<- factor(survey$continentNum, levels= c(1:5), labels = c("Europe", "Africa", "Americas", "Oceania", "Asia")) 
survey$gender<- factor(survey$gender, levels=c(1,2), labels=c("Man", "Woman"))
class(survey$continentNum)
class(survey$gender)
class(survey$wtpTot)

#with ggplot
#survey$continentNum<- factor(survey$continentNum, levels= c(1:5), labels = c("Europe", "Africa", "Americas", "Oceania", "Asia")) 
#class(survey$continentNum)
#boxplot
ggplot(survey, aes(x = continentNum, y =wtpTot)) + geom_boxplot() + facet_wrap(~gender)
#lineplot
survey %>% group_by(continentNum, gender) %>% summarise(wtpTot_groups = mean(wtpTot)) -> surveyP
surveyP %>% ggplot() + aes(x = continentNum, y = wtpTot_groups, color = gender) + geom_line(aes(group = gender)) + geom_point() + 
  labs(title = "Gender & Continent on wtpTot", x = "Continent", y="(Mean) Willingness to Pay", color= "Sex")


