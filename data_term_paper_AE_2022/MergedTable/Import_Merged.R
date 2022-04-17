  library(quantmod) 

  rm()

# Symbols
  symbols <- c("AAPL", "BAC"  , "CMCSA", "CVX"  , "C"     , "DIS"  , "GE"  , "HD"   , "IBM"  , "JNJ"   , "JPM"  , "KO"   , "MCD" , "MRK"    , "MSFT",
               "ORCL", "PEP"  , "PFE"  , "PG"   , "QCOM"  , "SLB"  , "VZ"  , "WFC"  , "XOM"  , "^AEX"  , "^ATX" , "^AXJO", "^BFX", "^BFX"   , "^DJA", 
               "^DJI", "^FTAS", "^FTMC", "^FTSE", "^GDAXI", "^GSPC", "^HSI", "^IBEX", "^IXIC", "^MDAXI", "^N100", "^N225", "^NDX", "^OMXIPI", "^RUT", "^SMSI")
# Create new environment
  data_env <- new.env()
# Load symbols into data_env
  getSymbols(Symbols = symbols,  from = "2010-01-01", env = data_env, srs = "yahoo")
# Extract the close column from each object and merge into one xts object
  close_data <- do.call(merge, eapply(data_env, Cl))
# View the head of close_data
  head(close_data)
# Save Data as CSV file
  write.zoo(close_data, file = "close_data.csv", sep = ",")
















