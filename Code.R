library(quantmod)
library(rugarch)
library(tseries)
library(forecast)

#Downloading the data
MCD_data <- getSymbols("MCD", auto.assign = FALSE, from = "2010-01-01") #Post-covid period ? 

#Extracting the adjusted close price
MCD_adj_close <- MCD_data[, 6]

#Calculating logarithmic returns
MCD_log_returns <- log(MCD_adj_close) - log(lag(MCD_adj_close))
sum(is.na(MCD_log_returns)) #Only one NA due to differencing
MCD_log_returns <- na.omit(MCD_log_returns) #Removing the NA
plot(MCD_log_returns)

#Stationarity
adf.test(MCD_log_returns) #Stationary

#Histogram
hist(MCD_log_returns, breaks = "FD", border = F, col = "grey", main = "MCD log-returns", xlab = "", ylab = "") #Fat tails, high kurtosis

#ACF and PACF
par(mfrow=c(1,2))
acf(MCD_log_returns, xlab = "", ylab = "", main = "ACF") #Significant lag 1, 6, 7...
pacf(MCD_log_returns, xlab = "", ylab = "", main = "PACF") #Significant lag 1, 6, 7...

#Ljung-Box test
for (lag_order in c(4, 8, 12)) {
  print(Box.test(MCD_log_returns, type = "Ljung-Box", lag = lag_order)) #Super significant at all lags
} 

#AR(1)
ar1 <- Arima(MCD_log_returns, order = c(1, 0, 0))
summary(ar1)

#ACF and PACF
par(mfrow=c(1,2))
acf(ar1$residuals, xlab = "", ylab = "", main = "ACF") #Significant lag 6, 7, 8...
pacf(ar1$residuals, xlab = "", ylab = "", main = "PACF") #Significant lag 6, 7, 8...

#AR(2)
ar2 <- Arima(MCD_log_returns, order = c(2, 0, 0))
summary(ar2)

#ACF and PACF
par(mfrow=c(1,2))
acf(ar2$residuals, xlab = "", ylab = "", main = "ACF") #Significant lag 6, 7, 8...
pacf(ar2$residuals, xlab = "", ylab = "", main = "PACF") #Significant lag 6, 7, 8...

#Auto.arima
auto.arima(MCD_log_returns, stationary = T, ic = "aic") #ARMA (4,5)

#ARMA(4,5)
arma45 <- Arima(MCD_log_returns, order = c(4, 0, 5))

#ACF and PACF
par(mfrow=c(1,2))
acf(arma45$residuals, xlab = "", ylab = "", main = "ACF") #Only distant lags significant
pacf(arma45$residuals, xlab = "", ylab = "", main = "PACF") #Only distant lags significant

#Ljung-Box test
for (lag_order in c(4, 8, 12)) {
  print(Box.test(arma45$residuals, type = "Ljung-Box", lag = lag_order)) #Cannot reject null
}

#GARCH(1,1)
garch_spec <- ugarchspec(mean.model = list(armaOrder = c(4, 5)), variance.model = list(garchOrder = c(1, 1)))
garch_fit <- ugarchfit(garch_spec, MCD_log_returns)
garch_fit #Everything super significant, Specification tests OK, may want to try TARCH since negative bias is significant

#Volatility
dev.off()
plot(sigma(garch_fit))

#iGARCH(1,1)
igarch_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)), variance.model = list(model = "iGARCH", garchOrder = c(1, 1)))
igarch_fit <- ugarchfit(igarch_spec, MCD_log_returns)
igarch_fit

#Volatility
plot(sigma(igarch_fit))
