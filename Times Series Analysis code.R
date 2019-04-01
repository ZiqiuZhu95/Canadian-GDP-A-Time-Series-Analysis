setwd('C:/Users/Ziqiu/OneDrive/Documents/University Docs/Stat443')
#setwd('C:/Users/zhux9990/Downloads')
CONS.data <- read.csv("CONS_Canada.csv")
GDP.data <- read.csv("GDP_CONS_Canada.csv")

library("tseries", lib.loc="~/R/win-library/3.5")
library("fGarch", lib.loc="~/R/win-library/3.5")
library("rugarch", lib.loc="~/R/win-library/3.5")

##################################Question 1######################################
wt2 <- as.numeric(CONS.data$Seasonally.Adjusted)
wt1 <- as.numeric(GDP.data$GDP)
time = as.character(GDP.data$ï..)
plot(wt1, type="l", xlab = "Time", ylab = "Raw Values",
     ylim = c(min(min(wt1), min(wt2)) - 500, max(max(wt1),max(wt2))+500), col="blue", lwd=1)
lines(wt2, col="red", lwd = 1)
legend(20, 1200000, legend = c("GDP","CONS"),lty = c(1,1), lwd=c(2.5,2.5), col=c("blue","red"), cex=0.7 )

##################################Question 2#######################################
Xt = log(wt1)
plot (Xt, type = "l", xlab = "Time", ylab = "Decomposed Values")
############### How do we describe the data series? Increasing? ##############
#### Just what we see from the graph ###

TSeq = 1:length(Xt)
TSmodel = lm(Xt~TSeq)
summary(TSmodel)
#### Summarize TS Model

dx = diff(Xt)
DSmodel = lm(dx~1)
summary(DSmodel)
#### Summarize DS Model

Yt_TS = TSmodel$residuals
Yt_DS = DSmodel$residuals
plot(Yt_TS,type="l", xlab="t", ylab= "Yt values")
plot(Yt_DS,type="l", xlab="t", ylab= "Yt values")


#######################################Question 3######################################
# Function for BIC
BIC <- function(res, k, N) {
  bic <- log(sum(res^2) / N)
  bic <- bic + log(N) * k / N
  bic
}

#obtaining the bic for TS model
N_TS <- length(Yt_TS)
bic.array.TS <- rep(NA, 9)

for(ii in 0:9)
{
  model.arima.TS <- arima(Yt_TS, order = c(ii, 0, 0), include.mean=FALSE)
  res.arima.TS <- model.arima.TS$residuals
  bic.array.TS[ii+1] <- BIC(res.arima.TS, ii, N_TS)
}
which(bic.array.TS == min(bic.array.TS)) ## Gives 3, but since it's ii+1, we want AR(2)

# Puttng BIC_TS into a table
BIC_table_TS <- data.frame(cbind(c(0:9),bic.array.TS))
rownames(BIC_table_TS) <- c()
colnames(BIC_table_TS) <- c("p","BIC")
round(BIC_table_TS[2], digits=3)

#Obtaining the bic for DS model
N_DS = length(Yt_DS)
bic.array.DS = rep(NA,10) #10 lags are enough
for (ii in 0:9)
{
  model.arima.DS <- arima(Yt_DS, order = c(ii, 0, 0),include.mean=FALSE)
  res.arima.DS <- model.arima.DS$residuals
  bic.array.DS[ii+1] <- BIC(res.arima.DS, ii, N_DS)
}
which(bic.array.DS == min(bic.array.DS)) ## Gives 2, but since it's ii+1 we want AR(1)
#But the question asks us to do an AR(2) since AR(1) is too easy

# Puttng BIC_DS into a table
BIC_table_DS <- data.frame(cbind(c(0:9),bic.array.DS))
rownames(BIC_table_DS) <- c()
colnames(BIC_table_DS) <- c("p","BIC")
round(BIC_table_DS[2], digits=3)

#AR(2) for a DS Approach
DSModel.AR2 = arima(Yt_DS, order = c(2, 0, 0),include.mean=FALSE)
summary(DSModel.AR2)

#AR(2) is best for TS approach
TSModel.AR2 = arima(Yt_TS, order = c(2,0,0),include.mean=FALSE)
summary(TSModel.AR2)
#We have Yt = 1.333Yt-1 - 0.338Yt-2
#### The variances are in the diagonal, the summary gives a covariance matrix of the coefficients####
#Then the variances are 0.0049 and 0.00496 for the coefficients
sigma2.hat.TS = sum(TSModel.AR2$residuals^2)/N_TS # or just TS.Model.AR2$sigma2
sigma.hat.TS = sqrt(sigma2.hat.TS)
#sigma is 0.0199998

phi_TS = TSModel.AR2$coef #gives us the phi values)
rhos_TS <- ARMAacf(phi_TS,lag.max=2) ## Extracting the theoretical acf values for k=0,..,6. 
for (i in 3:8)
{
  rhos_TS = c(rhos_TS, phi_TS[1]*rhos_TS[i] + phi_TS[2]*rhos_TS[i-1])
}
gamma0_TS <- (sigma2.hat.TS/(1-phi_TS[1]*rhos_TS[2]-phi_TS[2]&rhos_TS[3]))
#sqrt(gamma0_TS) =0.008311591

#Moving average weights found recursively for k=0,1,..,6
psis_TS <- rep(NA,9)
psis_TS[1] <- 1
psis_TS[2] <- phi_TS[1]
for (i in 2:8)
{
  psis_TS[i+1] = phi_TS[1] * psis_TS[i] + phi_TS[2] * psis_TS[i-1] 
}


####################################################QUESTION 3 END####################################################

sigma2.hat.DS = DSModel.AR2$sigma2
sigma.hat.DS = sqrt(sigma2.hat.DS)

phi_DS = DSModel.AR2$coef #gives us the phi values)

#Moving average weights found recursively for k=0,1,..,9
psis_DS <- rep(NA,9)
psis_DS[1] <- 1
psis_DS[2] <- phi_DS[1]
for (i in 2:8)
{
  psis_DS[i+1] = phi_DS[1] * psis_DS[i] + phi_DS[2] * psis_DS[i-1] 
}

######################Question 4 #############################
# Use result from the book?# Say how you came up with it #

# Estimate Et[Yt+k] for the TS models
# Why do we use the YT residuals to estimate a DS model.
# finding expected values of E_T[Y_{t+k}] based on estimated phi values
lastYt_TS <- Yt_TS[length(Yt_TS)]
Etk_TS <- rep(NA,9)
Etk_TS[1] <- lastYt_TS
Etk_TS[2] <- phi_TS[1] * Yt_TS[length(Yt_TS)] + phi_TS[2] * Yt_TS[length(Yt_TS)-1]
for(i in 3:9){
  Etk_TS[i] <- phi_TS[1]*Etk_TS[i-1] + phi_TS[2] * Etk_TS[i-2]  
}
#Etk1_TS = Etk_TS[-1]
#Etk_TS = Etk_TS[-length(Etk_TS)]
#Xtk_TS = Etk1_TS - Etk_TS + TSmodel$coefficients[2] ##Xt TS
Xtk_TS = rep(NA,9)
Xtk_TS[1]= lastYt_TS - Yt_TS[length(Yt_TS)-1] + TSmodel$coefficients[2]
for (i in 2:9){
  Xtk_TS[i] = Etk_TS[i] - Etk_TS[i-1] + TSmodel$coefficients[2]
}
#Estimate Var[Yt+k]
VarEtk_TS = rep(NA, 9)
VarEtk_TS[1] = 0
VarEtk_TS[2] = sigma2.hat.TS
VarEtk_TS[3] = sigma2.hat.TS * (1 + (psis_TS[2] - psis_TS[1])^2)
for (i in 4:9) {
  VarEtk_TS[i] = sigma2.hat.TS * (1 + sum((psis_TS[2:(i-1)] - psis_TS[1:(i-2)])^2))
}

XtkCIlow_Ts = rep(NA, 9)
XtkCIhigh_Ts = rep(NA, 9)
for (i in 1:9)
{
  XtkCIlow_Ts[i] = Xtk_TS[i] - 1.96*sqrt(VarEtk_TS[i])
  XtkCIhigh_Ts[i] = Xtk_TS[i] + 1.96*sqrt(VarEtk_TS[i])
}
yLim <- range(XtkCIlow_Ts, XtkCIhigh_Ts)
plot(c(0:8), Xtk_TS, xlab="k", ylab = "forecast", type= "l", ylim = yLim, xlim = c(0,9))
par(new = TRUE)
plot(c(0:8), XtkCIhigh_Ts, col="red", lwd = 1, ylim = yLim, ylab = "", xlab="", type= "l", xlim = c(0,9))
par(new = TRUE)
plot(c(0:8), XtkCIlow_Ts, col="blue", lwd = 1, ylim = yLim, ylab = "", xlab="", type= "l", xlim = c(0,9))
legend(20, 1200000, legend = c("Lower CI","Upper CI"),lty = c(1,1), lwd=c(2.5,2.5), col=c("blue","red"), cex=0.7 )

#Estimate E(Y_t+k) For DS model

lastYt_DS = Yt_DS[length(Yt_DS)]
Etk_DS = rep(NA,9)
Etk_DS[1] = lastYt_DS
Etk_DS[2] = phi_DS[1] * Yt_DS[length(Yt_DS)] + phi_DS[2] * Yt_DS[length(Yt_DS)-1]
for(i in 3:9){
  Etk_DS[i] = phi_DS[1] * Etk_DS[i-1] + phi_DS[2] * Etk_DS[i-2]
}
Xtk_DS = Etk_DS + DSmodel$coefficients #### Xt DS estimate
VarEtk_DS = rep(NA,9)
VarEtk_DS[1] = 0
VarEtk_DS[2] = sigma2.hat.DS
for (i in 3:9){
  VarEtk_DS[i] = sigma2.hat.DS * sum(psis_DS[1:(i-1)]^2)}

EtkCIlow_DS = rep(NA, 9)
EtkCIhigh_DS = rep(NA, 9)
for (i in 1:9)
{
  EtkCIlow_DS[i] = Xtk_DS[i] - 1.96*sqrt(VarEtk_DS[i])
  EtkCIhigh_DS[i] = Xtk_DS[i] + 1.96*sqrt(VarEtk_DS[i])
}

#From TSmodel, Xt = 12.60655 + 0.00826t + Yt?
#From DSmodel, Xt = 0.00899 + Yt

DSyLim <- range(EtkCIlow_DS, EtkCIhigh_DS)
plot(c(0:8), Xtk_DS, xlab="K", ylab = "Forecast", type= "l", ylim = DSyLim, xlim = c(0,9))
par(new = TRUE)
plot(c(0:8), EtkCIhigh_DS, col="red", lwd = 1, ylim = DSyLim, ylab = "", xlab="", type= "l", xlim = c(0,9))
par(new = TRUE)
plot(c(0:8), EtkCIlow_DS, col="blue", lwd = 1, ylim = DSyLim, ylab = "", xlab="", type= "l", xlim = c(0,9))

######################################## Question 5 ##########################################
#Dickey Fuller Test
library("tseries", lib.loc="~/R/win-library/3.5")
adf.test(Xt, k =5)

# p-value is 0.2865 
# This means that it is difference stationary?
#Then we also have to analyze how to do a Dickey Fuller Test.

######################################## Question 6 ##########################################

########################################## QUESTION 6 DS#######################################################

#Should probably test for multiple ARMA(PQ) models from our Yt_TS because of the different cutoffs. Then we run multiple diagnostic tests
#Partial autocorrelation and autocorrelation functions plotted against k
acf(Yt_DS) #Doesn't really look damped
pacf(Yt_DS) # Not damped.. cut off is?
#We test for an AR(1) because it looks like that is a possible cutoff and we test for an MA(4) because that also looks like a possible cutoff
Q6DSmodel.AR1 = arima(Yt_DS, order = c(1,0,0), include.mean=FALSE)
DS.AR1res = Q6DSmodel.AR1$residuals
#We test for both diagnostics then

#Diagnostics for AR(1)
DS1.std.res = (DS.AR1res - mean(DS.AR1res))/sd(DS.AR1res)
hist(DS1.std.res, breaks=50, freq=F, xlab = "Standardized Residuals of the DS AR(1)")
curve(dnorm(x),col = "red",add=T)

N.DS.RES1 = length(DS1.std.res)
K3.DS1 = 1/(N.DS.RES1) * sum(DS1.std.res^3)
K4.DS1 = 1/(N.DS.RES1) * sum(DS1.std.res^4)
J.STAT.DS1 = N.DS.RES1 * ((K3.DS1^2/6 + (K4.DS1-3)^2/24))
chisq.crit = qchisq(0.95,2)
J.H0.TESTDS1 = (J.STAT.DS1 > chisq.crit)
J.H0.TESTDS1 #If True, Reject
#False

#box pierce test from A4
M_DS = ceiling(sqrt(length(DS.AR1res)))
Box.test(DS.AR1res, type = "Box-Pierce", lag=M_DS)
#returns - p-value

#Overfitting with r =4, we test for an AR(5)
Q6DSmodel.AR5 = arima(Yt_DS, order = c(5,0,0), include.mean=FALSE)
LR_DS = N_TS * log(Q6DSmodel.AR1$sigma2 / Q6DSmodel.AR5$sigma2) 
hyp.test.DS <- LR_DS > qchisq(p = 0.95, df = 2) # if true reject H_0
hyp.test.DS == TRUE

###Arch(6) test

N.DS1 <- length(DS.AR1res)
DS1.res.sq <- DS.AR1res^2
DS1.ARCH.model <- lm(DS1.res.sq[-(1:6)]~DS1.res.sq[-c((1:5), N.DS1)] 
                    +  DS1.res.sq[-c(1:4,(N.DS1-1),N.DS1)] + DS1.res.sq[-c(1:3, (N.DS1-2):N.DS1)] 
                    + DS1.res.sq[-c(1:2, (N.DS1-3):N.DS1)] + DS1.res.sq[-c(1, (N.DS1-4):N.DS1)] 
                    + DS1.res.sq[-((N.DS1-5):N.DS1)] - 1)
R.squared.DS1 <- summary(DS1.ARCH.model)$r.squared
DS1.ARCH.test.stat <- N.DS1 * R.squared.DS1
DS1.ARCH.test.crit <- qchisq(0.95,6)
DS1.ARCH.Null.Hypothesis <- (DS1.ARCH.test.stat > DS1.ARCH.test.crit) # if true reject H_0
DS1.ARCH.Null.Hypothesis
#True

#Diagnostics for the MA(4)
Q6DSmodel.MA4 = arima(Yt_DS, order=c(0,0,4), include.mean=FALSE)
DS.MA4res = Q6DSmodel.MA4$residuals

#standardized normals
DS2.std.res = (DS.MA4res - mean(DS.MA4res))/sd(DS.MA4res)
hist(DS2.std.res, breaks=50, freq=F, xlab = "Standardized Residuals of the DS AR(1)")
curve(dnorm(x),col = "red",add=T)

N.DS.RES2 = length(DS2.std.res)
K3.DS2 = 1/(N.DS.RES2) * sum(DS2.std.res^3)
K4.DS2 = 1/(N.DS.RES2) * sum(DS2.std.res^4)
J.STAT.DS2 = N.DS.RES2 * ((K3.DS2^2/6 + (K4.DS2-3)^2/24))
chisq.crit = qchisq(0.95,2)
J.H0.TESTDS1 = (J.STAT.DS2 > chisq.crit)
J.H0.TESTDS1 #If True, Reject

#box pierce test from A4
M_DS2 = ceiling(sqrt(length(DS.MA4res)))
Box.test(DS.MA4res, type = "Box-Pierce", lag=M_DS2)
#returns - p-value
#p-value 0.1197 do not reject

#Overfitting with r =4, we test for an MA(8)
Q6DSmodel.MA8 = arima(Yt_DS, order = c(0,0,8), include.mean=FALSE)
LR_DS = N_TS * log(Q6DSmodel.MA4$sigma2 / Q6DSmodel.MA8$sigma2) 
hyp.test.DS2 <- LR_DS > qchisq(p = 0.95, df = 2) # if true reject H_0
hyp.test.DS2 == TRUE
#TRUE

N.DS2 <- length(DS.MA4res)
DS2.res.sq <- DS.MA4res^2
DS2.ARCH.model <- lm(DS2.res.sq[-(1:6)]~DS2.res.sq[-c((1:5), N.DS2)] 
                    + DS2.res.sq[-c(1:4,(N.DS2-1),N.DS2)] + DS2.res.sq[-c(1:3, (N.DS2-2):N.DS2)] 
                    + DS2.res.sq[-c(1:2, (N.DS2-3):N.DS2)] + DS2.res.sq[-c(1, (N.DS2-4):N.DS2)] 
                    + DS2.res.sq[-((N.DS2-5):N.DS2)] - 1)
R.squared.DS2 <- summary(DS2.ARCH.model)$r.squared
DS2.ARCH.test.stat <- N.DS2 * R.squared.DS2
DS2.ARCH.test.crit <- qchisq(0.95,6)
DS2.ARCH.Null.Hypothesis <- (DS2.ARCH.test.stat > DS2.ARCH.test.crit) # if true reject H_0
DS2.ARCH.Null.Hypothesis
#TRUE

############################################# Question 6 TS ###########################################################
#Partial autocorrelation and autocorrelation functions plotted against k
acf(Yt_TS) #Definitely damped. PACF is damped implies we get an AR(P)
pacf(Yt_TS) # Not damped, cut off is 1 so we choose AR(1)
# Then we have an AR(1)
Q6TSmodel.AR1 = arima(Yt_TS, order = c(1,0,0), include.mean=FALSE)
Q6TSmodel.res <- Q6TSmodel.AR1$residuals
TS.std.res <- (Q6TSmodel.res - mean(Q6TSmodel.res))/sd(Q6TSmodel.res)
hist(TS.std.res, breaks = 50, freq = F, xlab = "Standardized Residuals of the TS AR(1)")
curve(dnorm(x),col = "red", add = T)
# R is the additional lags that you're overfitting so we overfit 4 lags

#Jarque-Bern
N.TS.Res = length(TS.std.res)
K3.TS = 1/(N.TS.Res) * sum(TS.std.res^3)
K4.TS = 1/(N.TS.Res) * sum(TS.std.res^4)
J.stat.TS = N.TS.Res * ((K3.TS^2/6) + (K4.TS-3)^2/24)
chisq.crit=qchisq(0.95,2)
J.H0.Test = (J.stat.TS > chisq.crit)
J.H0.Test #if True, reject

#box pierce test from A4
M_TS = ceiling(sqrt(length(Q6TSmodel.res)))
Box.test(Q6TSmodel.res, type = "Box-Pierce", lag=M_TS)
#returns - p-value

#Overfitting with r =4, we test for an AR(5)
Q6TSmodel.AR5 = arima(Yt_TS, order = c(5,0,0), include.mean=FALSE)
LR_TS = N_TS * log(Q6TSmodel.AR1$sigma2 / Q6TSmodel.AR5$sigma2) 
hyp.test.TS <- LR_TS > qchisq(p = 0.95, df = 2) # if true reject H_0
hyp.test.TS == TRUE
#Reject H0

###Arch(6) test

N <- length(Q6TSmodel.res)
TS.res.sq <- Q6TSmodel.res^2
TS.ARCH.model <- lm(TS.res.sq[-(1:6)]~TS.res.sq[-c((1:5), N)] 
                 +  TS.res.sq[-c(1:4,(N-1),N)] + TS.res.sq[-c(1:3, (N-2):N)] 
                 + TS.res.sq[-c(1:2, (N-3):N)] + TS.res.sq[-c(1, (N-4):N)] 
                 + TS.res.sq[-((N-5):N)] - 1)
R.squared <- summary(TS.ARCH.model)$r.squared
TS.ARCH.test.stat <- N * R.squared
TS.ARCH.test.crit <- qchisq(0.95,6)
TS.ARCH.Null.Hypothesis <- (TS.ARCH.test.stat > TS.ARCH.test.crit) # if true reject H_0

################################### Question 6 END #########################################

###################################### QUESTION 7 BEGINS #################################################
SP.data = read.csv("SP_data.csv")
Pt = as.numeric(SP.data$P)
#Regress log(PT) over delta + log(Pt-1)
Q7model = lm(log(Pt[-1]) ~ 1 + log(Pt[-length(Pt)]))
Q7res = Q7model$residuals
Q7acf = rep(NA,10)
for (i in 1:10){
  Q7acf[i] = sum(Q7res[-c((length(Q7res) - (i-1)):length(Q7res))] * Q7res[-c(1:i)])/sum(Q7res^2)
}
#Test if residuals are normally distributed
Q7.std.res = (Q7res - mean(Q7res))/sd(Q7res)
hist(Q7.std.res, breaks = 50, freq = F, xlab = "Standardized Residuals of S&P", main = "")
curve(dnorm(x),col = "red", add = T)
Q7.K3 = mean(Q7.std.res^3)
Q7.K4 = mean(Q7.std.res^4)
Q7J.stat = length(Q7.std.res) * ((Q7.K3^2/6) + (Q7.K4-3)^2/24)
Q7J.H0.Test = (Q7J.stat > chisq.crit)
Q7J.H0.Test #if True, reject

model.GARCH = garchFit(formula= ~garch(1,1), data=Q7res)
model.GARCH

############# Question 8 ####################

lnPt = log(Pt)
egarch11.spec = ugarchspec(variance.model=list(model="eGARCH",
                                               garchOrder=c(1,1)),
                           mean.model=list(armaOrder=c(0,0)))
egarch12.spec = ugarchspec(variance.model=list(model="eGARCH",
                                               garchOrder=c(1,2)),
                           mean.model=list(armaOrder=c(0,0)))
egarch21.spec = ugarchspec(variance.model=list(model="eGARCH",
                                               garchOrder=c(2,1)),
                           mean.model=list(armaOrder=c(0,0)))
egarch22.spec = ugarchspec(variance.model=list(model="eGARCH",
                                               garchOrder=c(2,2)),
                           mean.model=list(armaOrder=c(0,0)))
egarch11.fit = ugarchfit(egarch11.spec, lnPt)
egarch12.fit = ugarchfit(egarch12.spec, lnPt)
egarch21.fit = ugarchfit(egarch21.spec, lnPt)
egarch22.fit = ugarchfit(egarch22.spec, lnPt)
q8res = egarch22.fit@ fit$residuals
q8res.std = (q8res - mean(q8res))/sd(q8res)
#Jarque-Bern
N.q8 = length(q8res.std)
K3.Q8 = 1/(N.q8) * sum(q8res.std^3)
K4.Q8 = 1/(N.q8) * sum(q8res.std^4)
J.stat.Q8 = N.q8 * ((K3.Q8^2/6) + (K4.Q8-3)^2/24)
chisq.crit=qchisq(0.95,2)
J.H0.Test.Q8 = (J.stat.Q8 > chisq.crit)
J.H0.Test.Q8 #if True, reject

#box pierce test from A4
M_Q8 = ceiling(sqrt(length(q8res)))
Box.test(q8res, type = "Box-Pierce", lag=M_Q8)
#returns - p-value

#Overfitting with r =4, we test for an AR(5)
Q6TSmodel.AR5 = arima(Yt_TS, order = c(5,0,0), include.mean=FALSE)
LR_TS = N_TS * log(Q6TSmodel.AR1$sigma2 / Q6TSmodel.AR5$sigma2) 
hyp.test.TS <- LR_TS > qchisq(p = 0.95, df = 2) # if true reject H_0
hyp.test.TS == TRUE
#Reject H0

###Arch(6) test

N <- length(q8res)
Q8.res.sq <- q8res^2
Q8.ARCH.model <- lm(Q8.res.sq[-(1:6)]~Q8.res.sq[-c((1:5), N)] 
                    +  Q8.res.sq[-c(1:4,(N-1),N)] + Q8.res.sq[-c(1:3, (N-2):N)] 
                    + Q8.res.sq[-c(1:2, (N-3):N)] + Q8.res.sq[-c(1, (N-4):N)] 
                    + Q8.res.sq[-((N-5):N)] - 1)
R.squared <- summary(TS.ARCH.model)$r.squared
Q8.ARCH.test.stat <- N * R.squared
Q8.ARCH.test.crit <- qchisq(0.95,6)
Q8.ARCH.Null.Hypothesis <- (Q8.ARCH.test.stat > Q8.ARCH.test.crit) # if true reject H_0

