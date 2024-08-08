###-------------------- Acetylene reduction rates responses to element concentrations and diazotrophic bacteria in Racomitrium lanuginosum (Hedw.)Brid.--------------------###

# Regression analyses of variables associated to moss N2-fixation. Dennis Escolástico-Ortiz. 2022

#install.packages("lme4")
# install.packages("lmerTest")
library(lme4)
packageVersion("lme4")
library(lmerTest)
library(ggplot2)
#install.packages("car")
require(car)
require(MASS)


##-------------------- 1. Explore data for modeling - Whole data set --------------------##

# Import data 

setwd("C:/Users/escol/OneDrive - Université Laval/6_Doctorat en Biologie/1_Thesis/Chapter3_Microbiome")
BNF_metals<- read.delim ("Correlation_analyses/BNF_Metal_content.txt",header=T)
head(BNF_metals)
summary(BNF_metals)
var(BNF_metals$BNF)


# Transform to factors some variables and check the summary of the data
str(BNF_metals)
BNF_metals$Habitat<-as.factor(BNF_metals$Habitat)
BNF_metals$Plot<-as.factor(BNF_metals$Plot)
str(BNF_metals)
head(BNF_metals)
summary(BNF_metals)

# Density of BNF (whole data)
theme_bw()

ggplot(BNF_metals, aes(BNF)) +
  geom_density()

# Density of BNF per habitat
ggplot(BNF_metals, aes(BNF)) +
  geom_density(aes(fill = Habitat), alpha=0.6) 

# Check distribution of Random effects
ggplot(BNF_metals, aes(x = BNF)) + geom_density() + facet_wrap("Plot")




# Explore relationships between the response and predictor variables
qplot(X31P, BNF, data = BNF_metals, 
      geom = c("point", "smooth"))

# Remove outliers X31P (elements)
BNF_metals<- BNF_metals[-c(58,60),]

qplot(X51V, BNF, data = BNF_metals, 
      geom = c("point", "smooth"))

qplot(X56Fe, BNF,  data = BNF_metals, 
      geom = c("point", "smooth"))

qplot(X98Mo, BNF,  data = BNF_metals, 
      geom = c("point", "smooth"))

qplot( X39K, BNF, data = BNF_metals, 
      geom = c("point", "smooth"))

# Re-scale the predictor variables for the model
BNF_metals[,4:8] <- scale(BNF_metals[,4:8],center=TRUE,scale=TRUE)


### Check which distribution fits the datta better ###
# Normal distribution
par(mfrow=c(1,2))
qqp(BNF_metals$BNF,"norm")

# Delete outliers and missing data
#Outliers
BNF_metals<- BNF_metals[-c(3,51,55), ]
BNF_metals
qqp(BNF_metals$BNF,"norm")

# Log-normal distribution
qqp(BNF_metals$BNF,"lnorm")

# Negative binomial
# Add a 1 to avoid zeros.
BNF_metals$BNF.t<- BNF_metals$BNF + 1
BNF_metals$BNF.t<- round(BNF_metals$BNF.t)
nbinom1 <- fitdistr(BNF_metals$BNF.t, "Negative Binomial")
qqp(BNF_metals$BNF.t, "nbinom", size = nbinom1$estimate[[1]], mu = nbinom1$estimate[[2]])

# Gamma distribution 
gamma1 <- fitdistr(BNF_metals$BNF, "gamma")
qqp(BNF_metals$BNF, "gamma", shape = gamma1$estimate[[1]], rate = gamma1$estimate[[2]])


##-------------------- 2. Linear Mixed Models - BNF --------------------##
#BNF_metals[,1] <- log(BNF_metals[,1])

#lmm0 <- lmer(BNF ~  X31P + X51V + X56Fe + X98Mo + Habitat + (1|Plot), data = BNF_metals,
#             REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
#summary(lmm0)

# Check residual distribution
#plot(fitted(lmm0), residuals(lmm0), xlab = "Fitted Values", ylab = "Residuals")
#abline(h = 0, lty = 2)
#lines(smooth.spline(fitted(lmm0), residuals(lmm0),df=2))
# The residuals of the linear mixed model are not normally distributed !!!
# Not a goof model

##-------------------- 3. Generalized Linear Mixed Models - BNF --------------------##

# Generalized linear model 
glm0 <- glm(BNF ~  X31P + X51V + X56Fe + X98Mo + Habitat , data = BNF_metals,
            family=gaussian(link = "log"))

summary(glm0)

# Check residual distribution
plot(fitted(glm0 ), residuals(glm0 ), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(glm0 ), residuals(glm0 ),df=2))

# Shapiro test for normality or residuals
shapiro.test(residuals(glm0))
studres(glm0)

# Penalized quasilikelihood (PQL) is a flexible technique that can deal with non-normal data, unbalanced design, and crossed random effects.
# However, it produces biased estimates if your response variable fits a discrete count distribution, 
# like Poisson or binomial, and the mean is less than 5 - or if your response variable is binary.

# Generalized Linear Mixed Models via PQL
PQL1<- glmmPQL(BNF ~ X31P + X51V + X56Fe + X98Mo + Habitat, ~ 1|Plot,
               family=gaussian(link = "log"),data=BNF_metals,verbose = FALSE)

summary(PQL1)

PQL2<- glmmPQL(BNF ~ X31P +X56Fe + Habitat, ~ 1|Habitat/Plot,
               family=gaussian(link = "log"),data=BNF_metals,verbose = FALSE)

summary(PQL2)

# Check residual distribution
plot(fitted(PQL1), residuals(PQL1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(PQL1), residuals(PQL1),df=2))

# Shapiro test for normality or residuals
shapiro.test(residuals(PQL1))


# Generalized Linear Mixed Models using Template Model Builder ###

# install.packages('TMB', type = 'source')
library(glmmTMB)

TMB1<- glmmTMB(BNF ~ X31P + X51V + X56Fe + X98Mo + Habitat ,
               data=BNF_metals,family=gaussian(link = "log"))
summary(TMB1)

# Check residual distribution
plot(fitted(TMB1), residuals(TMB1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(TMB1), residuals(TMB1),df = 2))

# Check normality
qqnorm(residuals(TMB1))
qqline(residuals(TMB1))

# Shapiro test for normality or residuals
shapiro.test(residuals(TMB1))
studres(TMB1)


TMB2<- glmmTMB(BNF ~ X31P + X51V + X56Fe + X98Mo + Habitat + (1|Plot),
               data=BNF_metals,family=gaussian(link = "log"))

summary(TMB2)

# Check residual distribution. Check if variance differs over predicted values.
plot(fitted(TMB2), residuals(TMB2), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(TMB2), residuals(TMB2), df = 2))

# Check normality
qqnorm(residuals(TMB2))
qqline(residuals(TMB2))

# Shapiro test for normality or residuals
shapiro.test(residuals(TMB2))

# Delete extreme values
BNF_metals<- BNF_metals[-c(54,53), ]

# Generalized Linear Mixed Models via Laplace approximation
# The Laplace approximation is a special case of a parameter estimation method called Gauss-Hermite quadrature (GHQ),
# with one iteration.

Laplace1 <- glmer(BNF ~ X31P + X56Fe + X51V + X98Mo + Habitat + (1|Plot),
                 data=BNF_metals, gaussian(link = "log"), glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e4)))
                 
summary(Laplace1)

# Check residual distribution
plot(fitted(Laplace1), residuals(Laplace1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Laplace1), residuals(Laplace1),df=2))

# Shapiro test for normality or residuals
shapiro.test(residuals(Laplace1))
studres(Laplace1)

# Model without non-significant effects
Laplace2 <- glmer(BNF ~ X31P + X56Fe +X51V + X98Mo + Habitat + (1|Plot),
                  data=BNF_metals, gaussian(link = "log"), glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))

summary(Laplace2)



# Check for singularity
tt <- getME(Laplace1,"theta")
ll <- getME(Laplace1,"lower")
min(tt[ll==0])
# Not a problem

# Double-checking gradient calculations
# Extract pre-computed information:
derivs1 <- Laplace1@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))

# One general problem is that large scaled gradients are often associated with small absolute gradients:
# we might decide that we're more interested in testing the (parallel) minimum of these two quantities:

max(pmin(abs(sc_grad1),abs(derivs1$gradient)))


# What if we redo the calculations with numDeriv?
dd <- update(Laplace1,devFunOnly=TRUE)
pars <- unlist(getME(Laplace1,c("theta","fixef")))
grad2 <- grad(dd,pars)
hess2 <- hessian(dd,pars)
sc_grad2 <- solve(hess2,grad2)
max(pmin(abs(sc_grad2),abs(grad2)))


# Restart from previous model
ss <- getME(Laplace1,c("theta","fixef"))
Laplace3 <- update(Laplace1,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

Laplace4 <- update(Laplace1,start=ss,control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
summary(Laplace4)

# Try a different optimizer
library("numDeriv")
library("RCurl") ## to source() from Github


#gm_all <- allFit(Laplace1)

# Gauss-Hermite quadrature (GHQ)

GHQ1 <- glmer(BNF ~ X31P +  X51V + X56Fe + X98Mo + Habitat + (1|Plot),
             data=BNF_metals, gaussian(link = "log"), nAGQ = 25)

summary(GHQ1)

GHQ2 <- glmer(BNF ~ X31P +  X51V + X98Mo + Habitat + (1|Habitat:Plot),
              data=BNF_metals, gaussian(link = "log"), nAGQ = 25)

summary(GHQ2)

# Check residual distribution
plot(fitted(GHQ1), residuals(GHQ1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(GHQ1), residuals(GHQ1),df = 2))

# Shapiro test for normality or residuals
shapiro.test(residuals(GHQ1))
studres(GHQ1)

# Check normality
qqnorm(residuals(GHQ1))
qqline(residuals(GHQ1))




# Error in updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ) : 
# nAGQ > 1 is only available for models with a single, scalar random-effects term


##-------------------- 4. Modeling using means per plot - Explore data --------------------##

setwd("C:/Users/escol/OneDrive - Université Laval/6_Doctorat en Biologie/1_Thesis/Chapter3_Microbiome")
Data_means<- read.delim ("Correlation_analyses/data_cor_means_top_20.txt",header=T)
Data_means
head(Data_means)
summary(Data_means)

# Transform to factors some variables and check the summary of the data
str(Data_means)
Data_means$Habitat<-as.factor(Data_means$Habitat)
Data_means$Plot<-as.factor(Data_means$Plot)
str(Data_means)
head(Data_means)
summary(Data_means)

ggplot(Data_means, aes(BNF, fill=Habitat)) +
  geom_density(alpha=.5) 

ggplot(Data_means, aes(BNF)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") 


ggplot(Data_means, aes(BNF,P, label=Plot)) +
  geom_point()+ geom_text(hjust=0, vjust=0)
# Seems negative correlation

# Outliers X31P (elements) - Plot 380
Data_means<- Data_means[-c(7),]

ggplot(Data_means, aes(BNF, V)) +
  geom_point()

# Seems negative correlation

ggplot(Data_means, aes(BNF, Fe)) +
  geom_point()
# Seems negative correlation

ggplot(Data_means, aes(BNF, Mo)) +
  geom_point()
# Seems negative correlation

ggplot(Data_means, aes(BNF, K)) +
  geom_point()
# Seems there is no correlation

# Check distribution of BNF

# Normal distribution
par(mfrow=c(1,2))
qqp(Data_means$BNF,"norm")

# Shapiro-Wilk test
shapiro.test(Data_means$BNF)


# Delete outliers and missing data
#Outliers
Data_means <- Data_means[-c(7), ]
#BNF_metals

# Log-normal distribution
qqp(Data_means$BNF,"lnorm")


# Log transformed BNF
Data_means$BNF.log <- log(Data_means$BNF)

ggplot(Data_means, aes(BNF.log)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=0.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") 


##-------------------- 5. Linear  Models - BNF means --------------------##

# Linear models. No random effects due to the averaging of values per plot.

lm1 <- lm(BNF ~  P + V + Fe + Mo + Habitat, data = Data_means)

summary(lm1)

glm1 <- glm(BNF ~  P + V + Fe + Mo + Habitat, data = Data_means,
            family=gaussian(link = "log"))

summary(glm1)

# Re-scale the predictor variables for the model
Data_means[,5:9] <- scale(Data_means[,5:9])

lm2 <- lm(BNF ~  P + V + Fe + Mo + Habitat, data = Data_means)

summary(lm2)

# Check residual distribution
plot(fitted(lm1), residuals(lm1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(lm1), residuals(lm1)))

# Shapiro test for normality or residuals
shapiro.test(residuals(lm1))
studres(lm1)

##-------------------- 6. Linear  Models - BNF  top 20 nifH OTU abundance (means) --------------------##


ggplot(Data_means, aes(BNF, log(Faith_nifH))) +
  geom_point()

lm3 <- lm(BNF ~ log(Shannon_nifH) + log(Simpson_nifH) + log(Faith_nifH) , data = Data_means)

summary(lm3)

# Transform relative abundaces using log10
Data_means[, 16:25] <-log(Data_means[,16:25]) +1 

lm4 <- lm(BNF ~ Rhodomicrobium + Dolichospermum + Nostoc + Trichormus + Iningainematapete + 
            Azorhizobium + Calothrix + Fischerella + Anabaena + Mastigocladus, data = Data_means)

summary(lm4)

lm5 <- lm(BNF ~ Rhodomicrobium + Dolichospermum + Nostoc + Trichormus + 
            Calothrix + Fischerella + Anabaena, data = Data_means)

summary(lm5)


# Check residual distribution
plot(fitted(lm4), residuals(lm4), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(lm4), residuals(lm4),df = 2))

my_residuals4<-residuals(lm4)
my_fit4<-fitted(lm4)
plot(my_fit4,my_residuals4)

# Shapiro test for normality or residuals
shapiro.test(residuals(lm4))
studres(lm4)


# Check normality
qqnorm(residuals(lm4))
qqline(residuals(lm4))

##-------------------- 7. Generalized Linear  Models - BNF  top 20 nifH OTU abundance (means) --------------------##

library(glmmTMB)
packageVersion("glmmTMB")
TMB3<- glmmTMB(BNF ~ Anabaena + Azorhizobium + Calothrix + Dolichospermum + Fischerella + Iningainematapete +
                 Mastigocladus + Nostoc + Rhodomicrobium + Trichormus,
               data = Data_means,family=gaussian(link = "log"))
summary(TMB3)

# Check residual distribution. Check if variance differs over predicted values.
plot(fitted(TMB3), residuals(TMB3), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(TMB3), residuals(TMB3), df = 2))

# Check normality
qqnorm(residuals(TMB3))
qqline(residuals(TMB3))

# Shapiro test for normality or residuals
shapiro.test(residuals(TMB3))
