library(quantmod)
library(rugarch)
library(aTSA)
library(tseries)
library(ggpubr)
library(forecast)
library(lmtest)
library(moments)
library(ggplot2)
library(MTS)

#Downloading the data
MCD_data <- getSymbols("MCD", auto.assign = FALSE, from = "2010-01-01", to = "2022-05-01") #Starting 2010-01-01 as the provided data set but extending to April 2022

#Extracting the close price
MCD_close <- MCD_data[, 4]
rm(MCD_data) #Removing an unneccessary object from the environment

########################
### Data description ###
########################

#Close price summary statistics
length(MCD_close)
mean(MCD_close)
sd(MCD_close)
skewness(MCD_close)
kurtosis(MCD_close)

#Testing stationarity
adf.test(MCD_close)

#Calculating logarithmic returns
MCD_log_returns <- log(MCD_close) - log(lag(MCD_close))
sum(is.na(MCD_log_returns)) #Only one NA due to differencing
MCD_log_returns <- na.omit(MCD_log_returns) #Removing the NA

#Log-returns summary statistics
length(MCD_log_returns)
mean(MCD_log_returns)
sd(MCD_log_returns)
skewness(MCD_log_returns)
kurtosis(MCD_log_returns)

#Average annualized volatility
sqrt(250)*sd(MCD_log_returns)

#Normality test
jarque.bera.test(MCD_log_returns)

#Plotting close price and log returns
par(mfrow = c(2, 1))
plot(MCD_close, main = "MCD closing stock price", grid.col = NA)
#addEventLines(xts("Covid-19", as.Date("2020-01-20")), col = "blue", lwd = 2, pos = 1, offset = 2) #Vertical line in time of Covid-19
#hist(MCD_close, breaks = "FD", border = F, col = "grey", main = "MCD closing price", xlab = "", ylab = "")
plot(MCD_log_returns, main = "MCD logarithmic returns", grid.col = NA)
#addEventLines(xts("Withdrawal", as.Date("2022-03-08")), col = "blue", lwd = 2, pos = 2) #Vertical line in time of withdrawal

#Zooming in
dev.off()
plot(MCD_log_returns["2021/"])
addEventLines(xts("Invasion", as.Date("2022-02-24")), col = "red", lwd = 2, pos = 2) #Vertical line in time of invasion
addEventLines(xts("Withdrawal", as.Date("2022-03-08")), col = "blue", lwd = 2, pos = 4) #Vertical line in time of withdrawal

#Histogram of log-returns
hist(MCD_log_returns, breaks = "FD", border = F, col = "black", xlab = "", ylab = "", main = "") #Fat tails, high kurtosis

#Ploting variance
plot(MCD_log_returns^2)
addEventLines(xts("Withdrawal", as.Date("2022-03-08")), col = "blue", lwd = 2, pos = 2) #Vertical line in time of withdrawal
addEventLines(xts("Covid-19", as.Date("2020-01-20")), col = "red", lwd = 2, pos = 2, offset = 2)

#Stationarity
adf.test(MCD_log_returns) #Stationary

#ACF and PACF
par(mfrow=c(1,2))
#acf(MCD_log_returns, xlab = "", ylab = "", main = "ACF") #Significant lag 1, 6, 7...
#pacf(MCD_log_returns, xlab = "", ylab = "", main = "PACF") #Significant lag 1, 6, 7...
acf_ret <- ggAcf(MCD_log_returns)+
  theme_minimal()+
  ggtitle("ACF")+
  ylab("")+
  xlab("")
pacf_ret <- ggPacf(MCD_log_returns)+
  theme_minimal()+
  ggtitle("PACF")+
  ylab("")+
  xlab("")
ggarrange(acf_ret, pacf_ret, ncol = 1, nrow = 2)

#Ljung-Box test
for (lag_order in c(4, 8, 12)) {
  print(Box.test(MCD_log_returns, type = "Ljung-Box", lag = lag_order)) #Super significant at all lags
} 

######################
### ARMA modelling ###
######################

#Defining a function that estimates given ARMA, returns p-values of coefficients, AIC, and performs the Ljung-Box test. If needed, also plots ACF and PACF
arma_complete <- function(ar_oder, ma_order, plot_ACF_PACF = F) {
  arima_model <- Arima(MCD_log_returns, order = c(ar_oder, 0, ma_order)) #Estimate the model
  print("P-values:", quote = F)
  print(arima_model$coef/sqrt(diag(arima_model$var.coef))) #P-values of coefficients
  print(paste0("AIC: ", arima_model$aic), quote = F)
  for (lag_order in c(4, 8, 12)) {
    print(Box.test(arima_model$residuals, type = "Ljung-Box", lag = lag_order)) #Ljung-Box test
  }
  if (plot_ACF_PACF == T) {
    par(mfrow=c(1,2))
    acf(arima_model$residuals, xlab = "", ylab = "", main = "ACF") #Plot ACF
    pacf(arima_model$residuals, xlab = "", ylab = "", main = "PACF") #Plot PACF
  }
}

#AR(1)
arma_complete(1, 0, T) #AR(1) significant but many dependencies left

#AR(2)
arma_complete(2, 0) #AR(2) insignificant, dependencies in 8th and 12th lag

#ARMA(1,1)
arma_complete(1, 1) #MA(1) insignificant, dependencies in 8th and 12th lag

#ARMA(2,2)
arma_complete(2, 2) #AR(1) and AR(2) insignificant, both MAs significant, still dependencies left

#ARMA(0,3)
arma_complete(0, 3) #MA(3) insignificant, dependencies in 8th and 12th lag

#Inspecting the suggestions from auto.arima
auto.arima(MCD_log_returns, stationary = T, ic = "aic") #ARMA (3,2)
auto.arima(MCD_log_returns, stationary = T, ic = "bic") #ARMA (1,0)

#ARMA(3,2)
arma_complete(3, 2, T) #Only AR(3) insignificant, dependencies in 8th and 12th lag

#ARMA(4,3)
arma_complete(4, 3) #Everything significant but dependencies left

#ARMA(4,4)
arma_complete(4, 4, T) #Couple of terms insignificant, starting to remove dependencies

#ARMA(5,4)
arma_complete(5, 4, T) #Satisfactory results but the model is pretty complicated

#ARMA(4,5)
arma_complete(4, 5, T) #Satisfactory results and higher AIC than previous but still pretty complicated

#Comparing the models statistically
unrestr <- Arima(MCD_log_returns, order = c(4, 0, 5))
restr <- Arima(MCD_log_returns, order = c(1, 0, 0))
lrtest(unrestr, restr) #Reject the null of restricted model being better
hist(restr$residuals, breaks = "FD", border = F, col = "black")
jarque.bera.test(restr$residuals) #Reject null

#######################
### GARCH modelling ###
#######################

#Testing for arch-effects
arma45 <- arima(MCD_log_returns, order = c(4, 0 ,5))
arch.test(arma45) #Significant arch effects
Box.test(arma45$residuals^2, type = "Ljung-Box", lag = 8) #Confirmed

#Defining a function that estimates a specified GARCH model
garch_complete <- function(ar_order, ma_order, arch_order, garch_order) {
  garch_spec <- ugarchspec(mean.model = list(armaOrder = c(ar_order, ma_order)), variance.model = list(garchOrder = c(arch_order, garch_order)))
  garch_fit <- ugarchfit(garch_spec, MCD_log_returns)
  print(garch_fit) #Print the results
  return(garch_fit) #Store the results
}

#GARCH(1,1) with ARMA(4,5)
garch11arma45 <- garch_complete(4, 5, 1, 1) #AR(1), MA(1), AR(3) insignificant, rest super significant, specification tests ok

#GARCH(1,1) with ARMA(3,2) Trying more parsimonious model
garch11arma32 <- garch_complete(3, 2, 1, 1) #All ARMA terms insignificant

#GARCH(1,1) with ARMA(3,3) Another shot at a parsimonious model
garch11arma33 <- garch_complete(3, 3, 1, 1) #Everything significant, spec tests ok

#GARCH(1,1) with ARMA(1,0) Super parsimonious model
garch11arma10 <- garch_complete(1, 0, 1, 1) #Specification tests OK :D Negative sign bias significant => may wanna try TARCH

#Inspecting standardized residuals
#GARCH(1,1) with ARMA(4,5)
resid_norm_garch11arma45 <- residuals(garch11arma45, standardize = TRUE)
par(mfrow=c(1,2))
acf(resid_norm_garch11arma45, xlab = "", ylab = "", main = "ACF") #Almost no dependencies
pacf(resid_norm_garch11arma45, xlab = "", ylab = "", main = "PACF") #Only one dependency in a very distant lag
#GARCH(1,1) with ARMA(1,0)
resid_norm_garch11arma10 <- residuals(garch11arma10, standardize = TRUE)
acf(resid_norm_garch11arma10, xlab = "", ylab = "", main = "ACF") #Very few dependencies
pacf(resid_norm_garch11arma10, xlab = "", ylab = "", main = "PACF") #Two dependencies

#Testing normality
jarque.bera.test(resid_norm_garch11arma45) #Reject null of normality
jarque.bera.test(resid_norm_garch11arma10) #Reject null of normality :(

#Volatility
dev.off()
plot(sigma(garch11arma10))
#lines(sigma(garch11arma45), col = "red") #The lines are basically the same
addEventLines(xts("Withdrawal", as.Date("2022-03-08")), col = "red", lwd = 2) #Vertical line in time of withdrawal

#Zooming in
dev.off()
plot(sigma(garch11arma10)["2021/"])
addEventLines(xts("Invasion", as.Date("2022-02-24")), col = "red", lwd = 2, pos = 2) #Vertical line in time of invasion
addEventLines(xts("Withdrawal", as.Date("2022-03-08")), col = "blue", lwd = 2, pos = 4) #Vertical line in time of withdrawal

#TARCH(1,1) with ARMA(1,0)
tarch11arma10_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)))
tarch11arma10 <- ugarchfit(spec = tarch11arma10_spec, data = MCD_log_returns)
tarch11arma10 #Sign bias significant

#Checking TARCH residuals
resid_norm_tarch11arma10 <- residuals(tarch11arma10, standardize = TRUE)
par(mfrow=c(1,2))
acf(resid_norm_tarch11arma10, xlab = "", ylab = "", main = "ACF") #Almost no dependencies
pacf(resid_norm_tarch11arma10, xlab = "", ylab = "", main = "PACF") #2 dependencies at distant lags
Box.test(resid_norm_tarch11arma10, type = "Ljung-Box", lag = 20) #Cannot reject the null of no autocorrelation
jarque.bera.test(resid_norm_tarch11arma10) #Reject the null of normality

#TARCH(1,1) with ARMA(4,5)
tarch11arma45_spec <- ugarchspec(mean.model = list(armaOrder = c(4, 5)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)))
tarch11arma45 <- ugarchfit(spec = tarch11arma45_spec, data = MCD_log_returns)
tarch11arma45 #Everything super significant, specification tests ok

#Checking TARCH residuals
resid_norm_tarch11arma45 <- residuals(tarch11arma45, standardize = TRUE)
par(mfrow=c(1,2))
acf(resid_norm_tarch11arma45, xlab = "", ylab = "", main = "ACF") #Almost no dependencies
pacf(resid_norm_tarch11arma45, xlab = "", ylab = "", main = "PACF") #1 dependency at a distant lag
Box.test(resid_norm_tarch11arma45, type = "Ljung-Box", lag = 20) #Cannot reject the null on no autocorrelation
jarque.bera.test(resid_norm_tarch11arma45) #Reject the null of normality

#May wanna consider using the Student's distribution given the shape of the log-returns and the repeatadly rejected hypothesis of residuals' normality

#GARCH(1, 1) with ARMA(1,0) and Student's distribution
garch11arma10stud_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(garchOrder = c(1, 1)), distribution.model = "std")
garch11arma10stud <- ugarchfit(garch11arma10stud_spec, MCD_log_returns)
garch11arma10stud #Everything super significant, specification tests ok, negative sign bias significant => TARCH

#Inspecting the residuals
resid_norm_garch11arma10stud <- residuals(garch11arma10stud, standardize = TRUE)
par(mfrow=c(1,2))
acf(resid_norm_garch11arma10stud, xlab = "", ylab = "", main = "ACF") #Almost no dependencies
pacf(resid_norm_garch11arma10stud, xlab = "", ylab = "", main = "PACF") #1 or 2 dependencies at distant lags
Box.test(resid_norm_garch11arma10stud, type = "Ljung-Box", lag = 20) #Cannot reject the null on no autocorrelation
jarque.bera.test(resid_norm_garch11arma10stud) #Reject the null of normality

#TARCH(1, 1) with ARMA(1,0) and Student's distribution
tarch11arma10stud_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)), distribution.model = "std")
tarch11arma10stud <- ugarchfit(tarch11arma10stud_spec, MCD_log_returns)
tarch11arma10stud #Everything super significant, specification tests ok

#Inspecting the residuals
resid_norm_tarch11arma10stud <- residuals(tarch11arma10stud, standardize = TRUE)
par(mfrow=c(1,2))
acf(resid_norm_tarch11arma10stud, xlab = "", ylab = "", main = "ACF") #Almost no dependencies
pacf(resid_norm_tarch11arma10stud, xlab = "", ylab = "", main = "PACF") #1 or 2 dependencies at distant lags
Box.test(resid_norm_tarch11arma10stud, type = "Ljung-Box", lag = 20) #Cannot reject the null on no autocorrelation
jarque.bera.test(resid_norm_tarch11arma10stud) #Reject the null of normality

#Plotting volatility
dev.off()
plot(sigma(tarch11arma10stud))
lines(sigma(garch11arma10), col = "red") #The lines are basically the same

#Checking whether alpha + beta is smaller than one (stability condition)
tarch11arma10stud@fit$coef[4]+tarch11arma10stud@fit$coef[5]
#Testing by comparing it with the iGARCH model (as restr vs unrestr)
igarch11arma10stud_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "iGARCH", garchOrder = c(1, 1)), distribution.model = "std")
igarch11arma10stud <- ugarchfit(igarch11arma10stud_spec, MCD_log_returns)
LR <- -2*(likelihood(garch11arma10stud)-likelihood(igarch11arma10stud))
pchisq(LR,1) #Reject the null that restriction holds => no iGARCH

#Adding dummy for the withdrawal
withdrawal_dummy <- matrix(ifelse(index(MCD_log_returns) >= as.Date("2022-03-08"), 1, 0), ncol = 1)
tarch11arma10stud_dummy_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), external.regressors = withdrawal_dummy), distribution.model = "std")
tarch11arma10stud_dummy <- ugarchfit(tarch11arma10stud_dummy_spec, MCD_log_returns)
tarch11arma10stud_dummy #External regressor not significant :/

#Trying the same for covid
covid_dummy <- matrix(ifelse(index(MCD_log_returns) >= as.Date("2020-01-20"), 1, 0), ncol = 1) #First Covid-19 case in USA
tarch11arma10stud_covid_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), external.regressors = covid_dummy), distribution.model = "std")
tarch11arma10stud_covid <- ugarchfit(tarch11arma10stud_covid_spec, MCD_log_returns)
tarch11arma10stud_covid #Covid super significant => Let's change the topic? :D

###################
### Final model ###
###################
ext_regs <- cbind(covid_dummy, withdrawal_dummy)
colnames(ext_regs) <- c("Covid", "Withdrawal")
final_model_spec <-ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), external.regressors = ext_regs), distribution.model = "std")
final_model <- ugarchfit(final_model_spec, MCD_log_returns)
final_model

#Inspecting the residuals
Box.test(final_model@fit$residuals, type = "Ljung-Box", lag = 8) #Reject null

#Inspecting the standardized residuals
stand_resid_final_model <- residuals(final_model, standardize = TRUE)
par(mfrow=c(1,2))
acf(stand_resid_final_model, xlab = "", ylab = "", main = "ACF") #Almost no dependencies
pacf(stand_resid_final_model, xlab = "", ylab = "", main = "PACF") #1 or 2 dependencies at distant lags
Box.test(stand_resid_final_model, type = "Ljung-Box", lag = 20) #Cannot reject the null of no autocorrelation
Box.test(stand_resid_final_model^2, type = "Ljung-Box", lag = 20) #Cannot reject the null of no autocorrelation
jarque.bera.test(stand_resid_final_model) #Reject the null of normality
archTest(stand_resid_final_model, 20)
hist(stand_resid_final_model, breaks = "FD", border = F, col = "black")

#ACF and PACF
acf_res <- ggAcf(stand_resid_final_model)+
  theme_minimal()+
  ggtitle("ACF")+
  ylab("")+
  xlab("")
pacf_res <- ggPacf(stand_resid_final_model)+
  theme_minimal()+
  ggtitle("PACF")+
  ylab("")+
  xlab("")
ggarrange(acf_res, pacf_res, ncol = 1, nrow = 2)

#Testing for integration
igarch_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "iGARCH", garchOrder = c(1, 1), external.regressors = ext_regs), distribution.model = "std")
igarch <- ugarchfit(igarch_spec, MCD_log_returns)
LR <- -2*(likelihood(final_model)-likelihood(igarch))
pchisq(LR, 1) #Reject the null that restriction holds => no iGARCH
lrtest(final_model, igarch)
alpha_beta_sum <- final_model@fit$coef[4]+final_model@fit$coef[5] #0.87
log(2)/log(0.87)

#Plotting estimated volatility
final_model_no_covid_spec <-ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)), distribution.model = "std")
final_model_no_covid <- ugarchfit(final_model_no_covid_spec, MCD_log_returns)
plot(MCD_log_returns, main = "MCD log-returns", grid.col = NA)
lines(sigma(final_model_no_covid), col = "blue")
lines(sigma(final_model), col = "red")
addEventLines(xts("Covid-19", as.Date("2020-01-20")), col = "green", lwd = 2, pos = 2, srt=90) #Vertical line in time of Covid-19
addLegend("topleft", on = 1, legend.names = c("Log-returns", "GARCH volatility (baseline model)", "GARCH volatility (no event dummies)"), col = c("black", "red", "blue"), lty = 1, bty = "n", lwd = c(2,1,1))

#Shorter period
MCD_log_returns_cut <- MCD_log_returns["2021-12/"]
length(MCD_log_returns_cut) #104 obs
plot(MCD_log_returns_cut)
addEventLines(xts("Withdrawal", as.Date("2022-03-08")), col = "blue", lwd = 2, pos = 2) #Vertical line in time of withdrawal
adf.test(MCD_log_returns_cut) #Stationary
auto.arima(MCD_log_returns_cut, stationary = T, ic = "aic") #ARMA(0,0)
auto.arima(MCD_log_returns_cut, stationary = T, ic = "bic") #ARMA(0,0)
acf(MCD_log_returns_cut) #No dependencies
pacf(MCD_log_returns_cut) #No depencencies
hist(MCD_log_returns_cut, breaks = "FD", border = F, col = "black")
short_term_model_spec <-ugarchspec(mean.model = list(armaOrder = c(0, 0)), variance.model = list(model = "sGARCH", garchOrder = c(1, 1), external.regressors = matrix(ext_regs[(nrow(ext_regs)-104+1):nrow(ext_regs), 2], ncol = 1)))
short_term_model <- ugarchfit(short_term_model_spec, MCD_log_returns_cut)
short_term_model
stand_resid_short_term_model <- residuals(short_term_model, standardize = TRUE)
Box.test(stand_resid_short_term_model, type = "Ljung-Box", lag = 10) #Cannot reject the null of no autocorrelation
Box.test(stand_resid_short_term_model^2, type = "Ljung-Box", lag = 10) #Cannot reject the null of no autocorrelation
jarque.bera.test(stand_resid_short_term_model) #Reject the null of normality
archTest(stand_resid_short_term_model, 10)

#########################
### Robustness checks ###
#########################

### ARMA(4,5) ###
robust_model_spec <-ugarchspec(mean.model = list(armaOrder = c(4, 5)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), external.regressors = ext_regs), distribution.model = "std")
robust_model <- ugarchfit(robust_model_spec, MCD_log_returns)
robust_model #Results basically the same

#Inspecting the standardized residuals
stand_resid_robust_model <- residuals(robust_model, standardize = TRUE)
par(mfrow=c(1,2))
acf(stand_resid_robust_model, xlab = "", ylab = "", main = "ACF") #Almost no dependencies
pacf(stand_resid_robust_model, xlab = "", ylab = "", main = "PACF") #1 dependency at a distant lag
Box.test(stand_resid_robust_model, type = "Ljung-Box", lag = 20) #Cannot reject the null of no autocorrelation
Box.test(stand_resid_robust_model^2, type = "Ljung-Box", lag = 20) #Cannot reject the null of no autocorrelation
jarque.bera.test(stand_resid_robust_model) #Reject the null of normality
archTest(stand_resid_robust_model, 20) #No ARCH effects
hist(stand_resid_robust_model, breaks = "FD", border = F, col = "black") #Looks good, bruh

#Comparing volatilities of final and robust
dev.off()
plot(sigma(final_model), lwd = 3, main = "Volatility", grid.col = NA)
lines(sigma(robust_model), col = "red") #Indistinguishable
addLegend("topleft", on = 1, legend.names = c("Final model volatility (TARCH(1,1) ARMA(1,0))", "Robustness check model volatility (TARCH(1,1) ARMA(4,5))"), col = c("black", "red"), lty = 1, bty = "n", lwd = c(3, 1))

### Normal distribution (rather than t)
robust_norm_model_spec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), external.regressors = ext_regs))
robust_norm_model <- ugarchfit(robust_norm_model_spec, MCD_log_returns)
robust_norm_model #Same as final