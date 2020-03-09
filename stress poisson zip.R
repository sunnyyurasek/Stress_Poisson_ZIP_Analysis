library(gridExtra)
library(dplyr)
library(ggplot2)
library(lessR)
library(MASS)
library(car) 

stress <- read.csv(file.path("C:/Users/syurasek/OneDrive - Constellation Brands/Documents/Northwestern/PREDICT 410/Computational 6", "STRESS.csv"),head=TRUE,sep=",")

summary(stress)
SummaryStats(stress)
head(stress)
mydata<-data.frame(stress)

Histogram(STRESS, data=mydata)
qqnorm(mydata$STRESS, pch = 1, frame = FALSE)
qqline(mydata$STRESS, col = "steelblue", lwd = 2)

###Creating 2 additional varibles since we have large number of zero in STRESS
mydata$stressed<-ifelse(mydata$STRESS>0,1,0) #stressed yes or no 
mydata$level<-ifelse(mydata$stressed==1, mydata$STRESS,NA) #if not stressed mark it as NA
head(mydata)

###OLS modeling
summary(olsmodel<-lm(formula = STRESS ~ COHES + ESTEEM + GRADES + SATTACH , data = mydata))


plot(formula = STRESS ~ COHES + ESTEEM + GRADES + SATTACH , data = mydata, main = "OLS stress model" )
abline(olsmodel, col='red')

layout(matrix(c(1,2,3,4),2,2)) 
plot(olsmodel)

outlierTest(olsmodel) 

###predict using OLS
mydata$y_hat1<-fitted(olsmodel)
Histogram(y_hat1)


ggplot(mydata, aes(sample = y_hat1)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Residual Distribution")

ggplot(mydata, aes(STRESS, y_hat1)) +
  geom_point() +
  stat_smooth(method = "lm") +
  labs( title = "Predicted Stress vs Actual Stress")

###log transformation has been done over value that's not 0 by adding 0.00001. 

mydata$logSTRESS <- log(mydata$STRESS+0.0001)

#mydata$logstress <- log(mydata$level) (or this can be done using level which is all data shows stress that's not 0)

summary(olsmodel2 <- lm(logSTRESS ~ COHES + ESTEEM + GRADES + SATTACH, data = mydata))

layout(matrix(c(1,2,3,4),2,2)) 
plot(olsmodel2)

mydata$y_hat2<-fitted(olsmodel2)
Histogram(y_hat2=NULL, data=mydata, bin.start=0,bin.width=2)
Histogram(y_hat2)

ggplot(mydata, aes(logSTRESS, y_hat2)) +
  geom_point() +
  stat_smooth(method = "lm") +
  labs( title = "Log Predicted Stress vs Actual Stress")

ggplot(mydata, aes(sample = y_hat2)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "OLS 2 Residual Distribution")


# try log OLS again by using level as Y
mydata$loglevel <- log(mydata$level)
summary(olsmodel2_1 <- lm(loglevel ~ COHES + ESTEEM + GRADES + SATTACH, data = mydata))

layout(matrix(c(1,2,3,4),2,2)) 
plot(olsmodel2_1)


###Poisson model
summary(poi1 <- glm(STRESS ~ COHES + ESTEEM + GRADES + SATTACH, data = mydata, family = "poisson"))
mydata$y_hatpoi<- fitted(poi1)
mydata$respoi <- residuals(poi1, type = "deviance")
Histogram(y_hatpoi, data=mydata)

#Poisson residual and predicted value evaluation 
ggplot(mydata, aes(respoi, fill = ..count..)) +
  geom_histogram()+
  labs(title = "Poisson Residual Histogram")

par(mfrow=c(2,2))
Plot(y_hatpoi, respoi, out_cut=.10)
Plot(y_hatpoi, respoi/sqrt(y_hatpoi))
Plot(respoi)
Plot(respoi/sqrt(y_hatpoi), out_cut=2, fences=TRUE, vbs_mean=TRUE)

Plot(y_hatpoi, STRESS)
Plot(log(y_hatpoi+1),log(STRESS+1))

###over-dispersed poisson model

summary(poi2 <- glm(STRESS ~ COHES + ESTEEM + GRADES + SATTACH, data = mydata, family = "quasipoisson"))
mydata$y_hatpoiq<- fitted(poi2)
mydata$respoiq <- residuals(poi2, type = "deviance")
Histogram(y_hatpoiq, data=mydata)

#Quasi-Poisson residual and predicted value evaluation 
ggplot(mydata, aes(respoiq, fill = ..count..)) +
  geom_histogram()+
  labs(title = "QuasiPoisson Residual Histogram")

###nbd model

summary(nbd1 <- glm.nb(STRESS ~ COHES + ESTEEM + GRADES + SATTACH, data = mydata, maxit = 100000))
mydata$y_hatnbd1<- fitted(nbd1)
mydata$respoinbd1 <- residuals(nbd1, type = "deviance")
Histogram(y_hatnbd1, data=mydata)

#NBD residual and predicted value evaluation 
ggplot(mydata, aes(respoinbd1, fill = ..count..)) +
  geom_histogram()+
  labs(title = "Negative Binomial Residual Histogram")

###Unobserved Heterogeneity using theta value from nbd1
v = 1/nbd1$theta
x = seq(0.0, 3, 0.1)
gd = data.frame(x, g = dgamma(x, shape = 1/v, scale = v))
ggplot(gd, aes(x, g)) + geom_line()
ggsave("gamdenr.png", width=500/72, height=400/72, dpi=72)

qgamma((1:3)/4, shape = 1/v, scale = v)

#quick model comparison between poisson and nbd
se <- function(model) sqrt(diag(vcov(model)))
round(data.frame(
  poisson.coef=coef(poi1),quasi.coef=coef(poi2),neg.binomial.coef=coef(nbd1),
  se.poisson=se(poi1),se.quasi=se(poi2),se.neg.binomial=se(nbd1)),4)

###AIC BIC comparison aginst three models

c(poisson_AIC=AIC(poi1), quasi_AIC=AIC(poi2), nbd_AIC=AIC(nbd1))

c(poisson_BIC=BIC(poi1), quasi_BIC=BIC(poi2), nbd_BIC=BIC(nbd1))

###family cohesion variable change

# Mean & SD

meanCOHES <- mean(mydata$COHES)
sdCOHES <- sd(mydata$COHES)
COHESlow<-meanCOHES-sdCOHES
COHESlow

COHEShigh<-meanCOHES+sdCOHES
COHEShigh

pred_COHES_low<-2.735 - 0.013*COHESlow- 0.024*0 - 0.23*0 - 0.016*0
pred_COHES_low

exp(pred_COHES_low)

pred_COHES_high<-2.735 - 0.013*COHEShigh- 0.024*0 - 0.23*0 - 0.016*0
pred_COHES_high

exp(pred_COHES_high)

#What is the expected percent difference in the number of stressful events for those at high and low levels of family cohesion?
(exp(pred_COHES_high)-exp(pred_COHES_low))/exp(pred_COHES_low)

# Calculate the groups
mydata <- mydata %>%
  mutate(groupCOHES = case_when(
    COHES < (meanCOHES - sdCOHES) ~ "Low",
    between(COHES, meanCOHES - sdCOHES, meanCOHES + sdCOHES) ~ "Middle",
    COHES > (meanCOHES + sdCOHES) ~ "High"
  ))

table(mydata$groupCOHES)

###logistic regression using stressed variable

summary(logit1 <- glm(stressed ~ COHES + ESTEEM + GRADES + SATTACH, family = binomial, data = mydata))

#interpretation by using odds direction transformation 
odds_dir1 <- round(exp(-0.020733) - 1, digits = 3)#COHES
odds_dir1

odds_dir2 <- round(exp(-0.018867) - 1, digits = 3)#ESTEEM
odds_dir2

odds_dir3 <- round(exp(-0.025492 ) - 1, digits = 3)#GRADES
odds_dir3

odds_dir4 <- round(exp(-0.027730 ) - 1, digits = 3)#SATTACH
odds_dir4

mydata$y_hatlogit<-fitted(logit1)

mydata$logit_y<- 3.516735-0.020733*mydata$COHES -0.018867 *mydata$ESTEEM - 0.025492 *mydata$GRADES- 0.027730 *mydata$SATTACH
mydata$pi<-exp(mydata$logit_y)/(1+exp(mydata$logit_y)) #probability of success
mydata$outcome<-ifelse(mydata$pi>0.5,1,0)

ggplot(mydata, aes(stressed, pi)) +
  geom_point() +
  labs(x = "X", y = "P(Success)", title = "Probability of Stress due to COHES-ESTEEM-GRADES-SATTACH")

Histogram(pi)

###ZIP model by hand adding AGE, and do step by step varaible evaluation 

summary(zip1 <- glm(stressed ~ COHES + ESTEEM + GRADES + SATTACH + AGE, family = binomial, data = mydata))

summary(zip2 <- glm(stressed ~ COHES, family = binomial, data = mydata))

summary(zip3 <- glm(stressed ~ COHES + ESTEEM, family = binomial, data = mydata))

summary(zip4 <- glm(stressed ~ COHES + ESTEEM + GRADES, family = binomial, data = mydata))

#zip 3 ranked #2 comparing to our previous logit1 model. almost identical AIC. so proceed further with zip3. 

mydata$y_hatzip3<- 3.06456-0.02539*mydata$COHES -0.03238 *mydata$ESTEEM
mydata$pi_zip3<-exp(mydata$y_hatzip3)/(1+exp(mydata$y_hatzip3)) #probability of success
Histogram(pi_zip3)

mydata$outcome_zip3<-ifelse(mydata$pi_zip3>0.5,1,0)

ggplot(mydata, aes(stressed, pi_zip3)) +
  geom_point() +
  labs(x = "X", y = "P(Success)", title = "Probability of Stress due to COHES-ESTEEM")

#continue model second part  by using level as the response varaible. compare between poisson and negative binomial distribution

summary(zip_poi<-glm(level ~ COHES + ESTEEM + GRADES + SATTACH, data = mydata, family = "poisson"))

summary(zip_poi1<-glm(level ~ COHES + ESTEEM + GRADES + SATTACH + AGE, data = mydata, family = "poisson"))

summary(zip_poi2<-glm(level ~ COHES, data = mydata, family = "poisson"))

summary(zip_poi3<-glm(level ~ COHES + ESTEEM, data = mydata, family = "poisson"))

summary(zip_poi4<-glm(level ~ COHES + ESTEEM + GRADES, data = mydata, family = "poisson")) #this is the best model with lowest AIC

summary(zip_nbd<-glm.nb(level ~ COHES + ESTEEM + GRADES, data = mydata, maxit = 100000))#did not beat poisson model in terms of AIC

mydata$y_hat_zip_poi4<- 2.254807-0.007263*mydata$COHES -0.021411 *mydata$ESTEEM-0.017447*mydata$GRADES
mydata$pred_zip_poi4<-exp(mydata$y_hat_zip_poi4)#poisson is the log funtion of Y

#combination model
mydata$y_hat_comb<-ifelse(mydata$pi_zip3<0.5, 0, mydata$pred_zip_poi4) #logistic prob success less than 50% then use poisson model.

Histogram(y_hat_comb)


###pscl model
install.packages("pscl")
library(pscl)

summary(zero1<-zeroinfl(formula = STRESS~ COHES + ESTEEM + GRADES + SATTACH | COHES + ESTEEM + GRADES + SATTACH, data=mydata))

summary(zero2<-zeroinfl(STRESS~ COHES + ESTEEM + GRADES + SATTACH | COHES + ESTEEM + GRADES + SATTACH, data = mydata, dist="negbin", EM = TRUE))

AICzero1<-AIC(zero1)
AICzero1

AICzero2<-AIC(zero2)
AICzero2


vuong(zero1, zero2)

summary(zero3<-zeroinfl(STRESS~ COHES + ESTEEM | COHES + ESTEEM + GRADES , data = mydata))

summary(zero4<-zeroinfl(STRESS~ COHES + ESTEEM | COHES + ESTEEM + GRADES , data = mydata, dist="negbin", EM = TRUE))

AICzero3<-AIC(zero3)
AICzero3

AICzero4<-AIC(zero4)
AICzero4

vuong(zero2, zero4)
