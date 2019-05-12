## Currency from Oanda
library(quantmod)
ts.index <- NULL

index.name <- "EUR/CHF"

ticker <- "EUR/CHF"

Sys.setenv(tz = "UTC")

from = "1998-01-01"

to = Sys.Date()

for (yy in substr(from, 1, 4) : substr(to, 1, 4))#y<-substr(from, 1, 4)
{

  tmp.from <- paste0(yy, "-01-01")

  tmp.to <- paste0(yy, "-12-31")

  tmp.ts.index <- try(getSymbols(ticker, src="oanda", from=tmp.from, to=tmp.to,

                                 freq = "daily", auto.assign = FALSE), silent=TRUE)

  if ("try-error" %in% class(tmp.ts.index)){

    next

  } else {

    print(paste0("Data downloaded for year: ", yy))

    ts.index <- rbind(ts.index, tmp.ts.index)

  }

}



index.levels <- ts.index[-1, 1]



# Returns in pips

index.returns <- diff(ts.index[, 1])

index.returns <- index.returns[-1, ]

