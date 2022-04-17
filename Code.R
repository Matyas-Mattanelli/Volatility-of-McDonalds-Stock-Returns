library(quantmod)
library(rugarch)
library(tseries)

#Downloading the data
MCD_data <- getSymbols("MCD", auto.assign = FALSE, from = "2021-01-01") #Post-covid period ? 

#Extracting the adjusted close price
MCD_adj_close <- MCD_data[, 6]

#Calculating logarithmic returns
MCD_log_returns <- log(MCD_adj_close) - log(lag(MCD_adj_close))
sum(is.na(MCD_log_returns)) #Only one NA due to differencing
MCD_log_returns <- na.omit(MCD_log_returns) #Removing the NA
plot(MCD_log_returns)

#Stationarity
adf.test(MCD_log_returns) #Stationary

#ACF and PACF
par(mfrow=c(1,2))
acf(MCD_log_returns, xlab = "", ylab = "", main = "ACF") #No dependencies
pacf(MCD_log_returns, xlab = "", ylab = "", main = "PACF") #Signficant lag 17, maybe due to the volatility?

#GARCH(1,1)
garch_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)), variance.model = list(garchOrder = c(1, 1)))
garch_fit <- ugarchfit(garch_spec, MCD_log_returns)
garch_fit #Beta is almost 1, may wanna try iGARCH

#Volatility
dev.off()
plot(sigma(garch_fit))

#iGARCH(1,1)
igarch_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)), variance.model = list(model = "iGARCH", garchOrder = c(1, 1)))
igarch_fit <- ugarchfit(igarch_spec, MCD_log_returns)
igarch_fit

#Volatility
plot(sigma(igarch_fit))
