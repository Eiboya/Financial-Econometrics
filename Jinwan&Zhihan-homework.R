# ========================= #
# Question 1, Summary Stats #
# ========================= #
library(tidyverse)
library(kableExtra)
library(patchwork)
library(moments)
library(anytime)
library(hrbrthemes)
library(psych)
library(bayesforecast)
library(ggplot2)
library(plotly)
library(PerformanceAnalytics)
library(zoo)
library(FitAR)
library(fGarch)
library(data.table) # To transpose dataframe

spy <- read.delim("F:/Download/daily_rets_day_ret_RV_BV_RJ_2.txt", header=FALSE)
spy$V1 <- anydate(spy$V1)
names(spy)[1] <- "date"
names(spy)[2] <- "DR"
names(spy)[3] <- "RV"
plot_rd <- spy %>%
  ggplot(aes(x=date, y=DR)) +
  geom_area(fill="blue", alpha=0.5) +
  geom_line(color="blue") +
  ylab("daily return") +
  theme_ipsum_es()
plot_rd <- ggplotly(plot_rd)
plot_rd
spy %>%
  summarise(Caption="Summary Stats for Daily Return",
            MeanDR=mean(spy$DR),
            VarDR=var(spy$DR),
            N=n(),
            SkeDR=moments::skewness(spy$DR),
            KurDR=moments::kurtosis(spy$DR)-3)
diff_dr <- spy$DR-mean(spy$DR)
n <- length(diff_dr)
diff_dr_cub <- diff_drˆ3
s <- sum(diff_dr_cub)
skewness_dr <- (var(spy$DR)ˆ(3/2)*n)ˆ(-1)*s
skewness_dr
diff_df_4 <- diff_drˆ4
k <- sum(diff_df_4)
kurtosis_dr <- (var(spy$DR)ˆ(2)*n)ˆ(-1)*k-3

kurtosis_dr
ggacf(spy$DR) + theme_bw()
spy$DRsq <- spy$DRˆ2
plot_rd_sq <- spy %>%
  ggplot(aes(x=date, y=DRsq)) +
  geom_area(fill="yellow", alpha=0.5) +
  geom_line(color="yellow") +
  ylab("daily return Squred") +
  theme_ipsum_es()
plot_rd_sq <- ggplotly(plot_rd_sq)
plot_rd_sq
ggacf(spy$DRsq)+theme_bw()
m30 <- apply.rolling(as.ts(spy$DR), 30, FUN="mean")
v30 <- apply.rolling(as.ts(spy$DR), 30, FUN="var")
s30 <- apply.rolling(as.ts(spy$DR), 30, FUN="skewness")
k30 <- apply.rolling(as.ts(spy$DR), 30, FUN="kurtosis")
par(mfrow=c(2,2))
plot.ts(m30)
plot.ts(v30)
plot.ts(s30)
plot.ts(k30)
m120 <- apply.rolling(as.ts(spy$DR), 120, FUN="mean")
v120 <- apply.rolling(as.ts(spy$DR), 120, FUN="var")
s120 <- apply.rolling(as.ts(spy$DR), 120, FUN="skewness")
k120 <- apply.rolling(as.ts(spy$DR), 120, FUN="kurtosis")
par(mfrow=c(2,2))
plot.ts(m120)
plot.ts(v120)
plot.ts(s120)
plot.ts(k120)
m365 <- apply.rolling(as.ts(spy$DR), 365, FUN="mean")
v365 <- apply.rolling(as.ts(spy$DR), 365, FUN="var")
s365 <- apply.rolling(as.ts(spy$DR), 365, FUN="skewness")
k365 <- apply.rolling(as.ts(spy$DR), 365, FUN="kurtosis")
par(mfrow=c(2,2))
plot.ts(m365)
plot.ts(v365)
plot.ts(s365)
plot.ts(k365)

# ======================= #
# Question 2, (a) and (b) #
# ======================= #

# ============= Method 1 of MLE ================ #
R_t2 <- numeric()
condvar2 <- numeric()
12
omega2 <- 0.01
alpha2 <- 0.05
beta2 <- 0.6
psi16 <- 2
condvar2[1] <- (omega2+alpha2*psi16ˆ2) / (1-alpha2-beta2)
observation_generation <- function(n){
  for (i in 1:n){
    R_t2[i] <- rnorm(1,0,condvar2[i])
    condvar2[i+1] <- omega2 + alpha2*(R_t2[i]-psi16)ˆ2 + beta2*condvar2[i]
  }
  return(R_t2[250:n])
}
# MLE calculations that takes 2 arguments, psi and y
myfunction2 <- function(psi,y){
  mu = psi[1]
  omega = exp(psi[2]) # exp so that omega>0
  alpha = exp(psi[3])
  beta = exp(psi[4]) # alpha + beta < 1
  theta = exp(psi[5])
  if (alpha + beta > 1){
    LT = -Inf # in order to put alpha + beta < 1
  } else {
    n = length(y) # sample size
    ell = numeric(n) # vector of size n
    condvar = numeric(n) # vector containing sigma_tˆ2
    condvar[1] = (omega+alpha*thetaˆ2)/(1-alpha-beta)
    for (t in 2:n){
      ell[t] = dnorm(y[t], mean = 0, sd=sqrt(condvar[t-1]),log=T)
      condvar[t] = omega + alpha * (y[t-1]-theta)ˆ2 + beta * condvar[t-1]
    }
    LT = sum(ell)
  }
  return(-LT)
}
psi1 = c(0, log(0.01), log(0.05), log(0.6), log(2))
# Replications, 100 times first to try
RR <- data.frame(replicate(n=5, optim(par=psi1, myfunction2,
                                      y=observation_generation(1000),
                                      control=list(maxit=9000000))$par))
#n_seq <- 1:2
#colnames(RR) <- paste("n", n_seq, sep = "_")
#rownames(RR) <- c("mu", "omega", "alpha", "beta", "theta")
#RR <- transpose(RR)

# ============ Method 2 ============ #
Garch_DG <- function (n, mu, omega, alpha, beta, theta) {
  burn = 250
  n = n + burn
  z = d.g(n)
  h = rep(0, n)
  r = rep(0, n)
  for (t in 2:n) {
    h[t] = omega + alpha * r[t-1]ˆ2 + beta * h[t-1]
    r[t] = mu + z[t] * sqrt(h[t])
  }
  return(r[burn:n])
}
# Run the following to check the generation
# r1 <- Garch_DG(1000, 0, 0.01, 0.05, 0.9, 0)
garchInit = function(repli) {
  Mean = mean(repli)
  Var = var(repli)
  S = 1e-6
  psi1 = c(mu = Mean, omega = 0.1*Var, alpha = 0.1, beta = 0.9, theta = 0)
  lowerBounds = c(mu = -10*abs(Mean),
                  omega = Sˆ2, alpha = S, beta = S, theta = -3)
  upperBounds = c(mu = 10*abs(Mean),
                  omega = 100*Var, alpha = 1-S, beta = 1-S, theta = 3)
  cbind(psi1=psi1, lowerBounds=lowerBounds, upperBounds=upperBounds)
}
# garchInit(r1)
# garchInit(r1)[,1]
Garch_NLL = function(psi, repli){
  mu = psi[1]; omega = psi[2]; alpha = psi[3]; beta = psi[4]; theta = psi[5]
  z = (repli - mu)
  Mean = mean((zˆ2))
  e = omega + alpha*c(Mean, (z[-length(repli)] - theta)ˆ2)
  h = timeSeries::filter(e, beta, "r", init = Mean)
  hh = sqrt(abs(h))
  llh = -sum(log(g.d(z, hh)))
}
Garch_Optim = function(repli, psi1, lowerBounds, upperBounds){
  fit = nlminb(start = psi1, objective = Garch_NLL,
               lower = lowerBounds,
               upper = upperBounds,
               control = list(trace=10),
               repli = repli)
}
# Change below for different number of replications
R <- 10
# Change below for different number of observations
O <- 100
# Change below for different designs
# Design Indicator D, D=1 for first and otherwise the second
D <- 2
# Change below for different distribution of zt
# Indicator for dist. of zt (˜D(0,1)/˜T(10), therefore "norm" or "std")
dist = "norm"
# Distribution for later MLEs
Garch_Dist <- dist
# Change below for different choice of variance process
vp <- 2
# parameter subcases
if (D == 1){
  mu = 0
  omega = 0.01
  alpha = 0.05
  beta = 0.9
  theta = 0
}
if (D == 2){
  mu = 0
  omega = 0.01
  alpha = 0.05
  beta = 0.6
  theta = 2
}
# distribution of zt, subcases
if (Garch_Dist == "norm"){
  d.g = function(n){
    rnorm(n,0,1)
  }
}
if (Garch_Dist == "std"){
  d.g = function(n){
    rt(n, df=10)
  }
}
# variances process subcases
if (vp == 1){
  v_eq = function(omega, alpha, Mean, z, repli){
    omega + alpha*c(Mean, (z[-length(repli)])ˆ2)
  }
}
if (vp == 2){
  v_eq = function(omega, alpha, Mean, z, repli, theta){
    omega + alpha*c(Mean, (z[-length(repli)] - theta)ˆ2)
  }
}
if (Garch_Dist == "norm"){
  g.d = function(z, hh){
    dnorm(x = z/hh)/hh
  }
}
if (Garch_Dist == "std"){
  g.d = function(z, hh){
    dt(x=z/hh, df=10)/hh
  }
}
mat = matrix(0, nrow = R, ncol = 5)
for (i in 1:R){
  r = Garch_DG(n = n,
               mu = mu, omega = omega, alpha = alpha,
               beta = beta, theta = theta)
  r = ts(r)
  RR1 = garchInit(r)
  psi1 = RR1[,1]
  lowerBounds = RR1[,2]
  upperBounds = RR1[,3]
  fit = Garch_Optim(r, psi1, lowerBounds, upperBounds)
  mat[i,] = fit$par
}
mest = colMeans(mat)
sdst = colSds(mat)
mest
sdst

# =================================== #
# 2.c Fit GARCH #
# 2.c.1 Model Selection #
# 2.c.2 Leverage Effect #
# #
# Garch_Model_Selection Function: #
# #
# With ARMA(1,1) Assumption; #
# Input takes a choice of garch.model #

# among baseline, TGARCH and EGARCH #
# and garch.order (1:2, 1:2) #
# =================================== #

Garch_Model_Selection <- function(garch.model, garch.order){
  mean.model = c(1,1);
  dist.choice = "norm";
  R = spy$DR;
  out.sample = length(R)/5;
  if (garch.model == "GARCH" | garch.model == "TGARCH") {
    spec.selection <- ugarchspec(variance.model = list(model = "fGARCH",
                                                       garchOrder=garch.order,
                                                       submodel = "GARCH"),
                                 mean.model = list(include.mean=TRUE,
                                                   armaOrder=mean.model),
                                 distribution.model = dist.choice)
    fit.outcome <- ugarchfit(spec = spec.selection,
                             data = R, out.sample = out.sample)
  }
  else if (garch.model == "eGARCH"){
    spec.selection = ugarchspec(variance.model = list(model = garch.model,
                                                      garchOrder=garch.order),
                                mean.model = list(include.mean=TRUE,
                                                  armaOrder=mean.model),
                                distribution.model = dist.choice)
    fit.outcome = ugarchfit(spec = spec.selection,
                            data = R, out.sample = out.sample)
  }
  return(fit.outcome)
}
Garch_Model_Selection("GARCH", c(1,1))
Garch_Model_Selection("GARCH", c(1,2))
Garch_Model_Selection("GARCH", c(2,1))
Garch_Model_Selection("GARCH", c(2,2))

# Testing for Leverage Effect, start using TGARCH and eGARCH
# Based on the previous selection result we chose GARCH(2,2)
Garch_Model_Selection("TGARCH", c(1,1))
Garch_Model_Selection("TGARCH", c(1,2))
Garch_Model_Selection("TGARCH", c(2,1))
Garch_Model_Selection("TGARCH", c(2,2))

# For TGARCH choose order (2,2)
Garch_Model_Selection("eGARCH", c(1,1))
Garch_Model_Selection("eGARCH", c(1,2))
Garch_Model_Selection("eGARCH", c(2,1))
Garch_Model_Selection("eGARCH", c(2,2))

# For eGARCH choose order (2,1)
Garch_Model_Selection("TGARCH", c(2,2))
Garch_Model_Selection("eGARCH", c(2,1))

# ============================= #
# Out of the Sample Diagnostics #
# ============================= #
model.t <- Garch_Model_Selection("TGARCH",
                                 c(2, 2))
model.e <- Garch_Model_Selection("eGARCH",
                                 c(2, 1))
model.t1 <- ugarchforecast(model.t,
                           data=NULL,
                           n.head=1,
                           n.roll=499,
                           out.sample=500)
model.t2 <- ugarchforecast(model.e,
                           data=NULL,
                           n.head=1,
                           n.roll=499,
                           out.sample=500)
length(spy$DR)/5
y = spy$DR[1998:length(spy$DR)]

yhat1 = model.t1@forecast$seriesFor[1,]
epsilon1 = y - yhat1
sigma.hat1 = model.t1@forecast$sigmaFor[1,]
z1 = epsilon1 /sigma.hat1
yhat2 = model.t2@forecast$seriesFor[1,]
epsilon2 = y - yhat2
sigma.hat2 = model.t2@forecast$sigmaFor[1,]
z2 = epsilon2 /sigma.hat2

Box.test(z1, lag=2, type = "Ljung-Box")
Box.test(z1, lag=6, type = "Ljung-Box")
Box.test(z1, lag=10, type = "Ljung-Box")
Box.test(z2, lag=2, type = "Ljung-Box")
Box.test(z2, lag=6, type = "Ljung-Box")
Box.test(z2, lag=10, type = "Ljung-Box")

# ============ #
# Risk Metrics #
# ============ #

rm.spec = ugarchspec(mean.model=list(armaOrder=c(1,1), include.mean=TRUE),
                     variance.model=list(model="TGARCH"),
                     distribution.model = "ged",
                     fixed.pars = list(omega = 0))
rm.model = ugarchfit(spec = rm.spec, data = spy$DR, out.sample = 500)
rm.pred = ugarchforecast(rm.model,
                         data=NULL,
                         n.head=1,
                         n.roll = 499,
                         out.sample = 500)

rm.pred = as.ts(rm.pred@forecast$sigmaFor[1,])ˆ2
rv.comp = as.ts(dat$RV[1998:2497])
plot(rm.pred, col="1", lty= 2, lwd = 1.5, main="Risk Metrics Versus Realised Volatility")
lines(rv.comp, col="2", lty=1.5, lwd = 1)
l1 <- as.expression("Risk Metrics")
l2 <- as.expression("Realised Volatility")
legend("topright",
       c(l1, l2),
       col=c(1, 2),
       lty=c(2, 1.5))

# ================================== #
# Question 3: Best AR and Comparison #
# ================================== #
library(FitAR)
library(forecast)
selection <- SelectModel(as.ts(spy$RV[1:1997]), lag.max = 14, ARModel = "AR",
                         Criterion = "AIC", Best = 6); selection
fit.1 <- arima(spy$RV, order = c(12,0,0))
fit.2 <- Arima(as.ts(spy$RV), model = fit1)
fc <- window(fitted(fit.2), start = 1998)
plot(as.ts(spy$RV), col = 4, main = "Model Choice AR(12)",
     ylab = "var", lty = 1)
lines(fc, col = 2, lty = 1)
legend("topright", c("Actual", "Prediction"), col=c(4,2),
       lty=c(1,2))

# ========================= #
# Question 4, Value-at-risk #
# ========================= #

#1. Three models and VaR when h=2. Similar way as that on moodle.

library(fGarch)
myfit = garchFit(˜garch(1,1), data=R, leverage=NULL, cond.dist="ged")
#myfit = garchFit(˜garch(1,1), data=R, leverage=TRUE, cond.dist="ged")
#myfit = garchFit(˜garch(1,1), data=R, cond.dist="norm")
results = myfit@fit
results = results$coef
cc = results[1]
omega = results[2]
alpha = results[3]
beta = results[4]
condvol = myfit@sigma.t # conditional volatility (sqrt(variance))
plot.ts(condvol)
# Simulation
S= 100000
yt = R[n]
sigmat = sqrt(omega + alpha * (yt-cc)ˆ2 + beta * condvol[n]ˆ2) # cond.vol for y_{t+1}
yforecast = numeric(S)
Xforecast = numeric(S)
set.seed(1230910)
for (s in 1:S){
  z = rnorm(n=2)
  y1= cc + sigmat * z[1]
  sigmatplusone = sqrt(omega + alpha * (sigmat * z[1])ˆ2 + beta * sigmatˆ2)
  y2 = cc + sigmatplusone * z[2]
  yforecast[s] = y1 + y2
  Xforecast[s] = 100 * 10 * (exp(y1+y2)-1)
}
hist(yforecast,20)
hist(Xforecast,30)
# Compute VaR of level 5%
q005 = quantile(yforecast, 0.05)
VaR_y = - q005
q005X = quantile(Xforecast, 0.05)
VaR_X = -q005X
# Compute ES of level 5%
subsety = yforecast[yforecast<q005]
length(subsety) #50 because 1000 * 0.05 = 50
hist(subsety,20)
ES_y = - mean(subsety)
subsetX = Xforecast[Xforecast < q005X]
length(subsetX) #5000 because 100000 * 0.05
hist(subsetX,20)
ES_X = - mean(subsetX)

# Results
c(VaR_y, ES_y)
c(VaR_X, ES_X)

# 2. VaR & Violations of the VaR

# Preparation of Functions and Specifications
spec1<-ugarchspec(variance.model = list(model = "fGARCH",
                                        garchOrder=c(2,2),
                                        submodel = "GARCH"),
                  distribution.model = "ged")
spec2<-ugarchspec(variance.model = list(model = "eGARCH",
                                        garchOrder=c(2,1),
                                        submodel = "GARCH"),
                  distribution.model = "ged")
spec3<-ugarchspec(variance.model = list(model = "fGARCH",
                                        garchOrder=c(2,2),
                                        submodel = "GARCH"),
                  distribution.model = "norm")
spec4<-ugarchspec(variance.model = list(model = "eGARCH",
                                        garchOrder=c(2,1),
                                        submodel = "GARCH"),
                  distribution.model = "norm")
spec5<-ugarchspec(variance.model = list(model = "fGARCH",
                                        garchOrder=c(2,2),
                                        submodel = "GARCH"),
                  distribution.model = "std")
spec6<-ugarchspec(variance.model = list(model = "eGARCH",
                                        garchOrder=c(2,1),
                                        submodel = "GARCH"),
                  distribution.model = "std")
hitseq <- function(myfit, alpha = 0.05) {
  VaR = -quantile(myfit, alpha)
  # Hit Sequence
  n = length(R[1:2497])
  hit = rep(0, n)
  for (i in 1:2497){
    if(R[i] < -VaR[i]){
      hit[i] = 1
    }
    
  }
  return(hit)
}

#Test 1: Uncondtional Ber(alpha) ####

test1uc = function(hit, alpha = 0.05){
  n = length(hit)
  T1 = length(which(hit==1))
  T0 = n - T1
  pi = T1 / n
  Lp = (1-pi)ˆT0 * piˆT1
  La = (1-alpha)ˆT0 * alphaˆT1
  #LR unconditional
  LRuc = 2*(log(Lp) - log(La))
  LRuc
}
#Test 2: Independence ####

test2ind = function(hit, alpha = 0.05){
  # T00
  n = length(hit)
  T00 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 0 & hit[i+1] == 0){
      T00 = T00 + 1
    }
  }
  # T01
  T01 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 0 & hit[i+1] == 1){
      T01 = T01 + 1
    }
  }
  # T10
  T10 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 1 & hit[i+1] == 0){
      T10 = T10 + 1
    }
  }
  # T11
  T11 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 1 & hit[i+1] == 1){
      T11 = T11 + 1
    }
  }
  # Likelihood of pi
  
  pi.01 = T01 / (T00 + T01)
  pi.11 = T11 / (T10 + T11)
  pi.00 = 1 - pi.01
  pi.10 = 1 - pi.11
  L.pi2 = (1 - pi.01)ˆT00 * pi.01ˆT01 * (1-pi.11)ˆT10 * pi.11ˆT11
  L.alpha2 = (1 - alpha)ˆT00 * alphaˆT01 * (1-alpha)ˆT10 * alphaˆT11
  LRcc = 2 * ( log(L.pi2) - log(L.alpha2) )
  LRcc
}

#Test 3: Add explanatory variables###

ttest <- function(reg, coefnum, val){
  tstat <- (coef(summary(reg))[coefnum,1]-val)/coef(summary(reg))[coefnum,2]
  pvalue <- 2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
  return(c(tstat, pvalue))
}
library(lmtest)
test3add = function(hit, set) {
  n = length(hit)
  T = length(data[,4])-1
  test3 = lm(hit[2:length(hit)] ˜ data[,4][1:2496] + data[,5][1:2496])
  print("b0=alpha")
  print(ttest(test3, 1, 0.05))
  print("b1 = 0")
  print(ttest(test3, 2, 0))
  print("b2 = 0")
  print(ttest(test3, 3, 0))
}
myfit1<-ugarchfit(spec=spec1, data = R)
myfit2<-ugarchfit(spec=spec2, data = R)
myfit3<-ugarchfit(spec=spec3, data = R)
myfit4<-ugarchfit(spec=spec4, data = R)
myfit5<-ugarchfit(spec=spec5, data = R)
myfit6<-ugarchfit(spec=spec6, data = R)

#Calculate VaR
VaR1 = -quantile(myfit1, 0.05)
VaR2 = -quantile(myfit2, 0.05)
VaR3 = -quantile(myfit3, 0.05)
VARTABLE<-cbind.data.frame(VaR1,VaR2,VaR3)
head(VARTABLE)
tail(VARTABLE)

#Hit Sequence
hit1<-hitseq(myfit1,0.05)
hit2<-hitseq(myfit2,0.05)
hit3<-hitseq(myfit3,0.05)
hit4<-hitseq(myfit4,0.05)
hit5<-hitseq(myfit5,0.05)
hit6<-hitseq(myfit6,0.05)

#Test 1
test1uc(hit1)
test1uc(hit2)
test1uc(hit3)
test1uc(hit4)
test1uc(hit5)
test1uc(hit6)

#Test 2
test2ind(hit1)
test2ind(hit2)
test2ind(hit3)
test2ind(hit4)
test2ind(hit5)
test2ind(hit6)
#Test 3
test3add(hit1)
test3add(hit2)
test3add(hit3)
test3add(hit4)
test3add(hit5)
test3add(hit6