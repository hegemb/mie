
##** MAKE DESCRIPTIVE TABLE ** 
#library(RGraphics) # support of the "R graphics" book, on CRAN
#library(gridExtra) 

## Make functions:
# Make function to table percent:
tabl1Perc <- function(lp,vNames,dd=data){ ## lp is lpnr, vNames is name of covariates.
  subdd <- dd[as.character(lp),] ## pick out the cases with lpnr = lp.
  sapply(subdd[,vNames],function(ss)round(prop.table(table(ss,useNA="no"))[2]*100,1)) 
  }
## Make function to table mean:
tabl1Mean <- function(x,vNames) apply(data[as.character(x),vNames],2,mean,na.rm=T)

## Make function to extract and create a table, using the two functions above:
descrTabl1 <- function(c.by.histo.LPNR, contVar=hc, dichoVar=op, allData=data, tabl1Mean.=tabl1Mean, tabl1Perc. = tabl1Perc, rNames){ #rNames is the names to be printed in the row of the table.
  t11.u <- round(sapply(c.by.histo.LPNR, function(x) tabl1Mean.(x,hc)),1)
  t12.u <- sapply(c.by.histo.LPNR, function(x)tabl1Perc.(lp = x,vNames=dichoVar))
  cohortPerc <- tabl1Perc.(lp = allData$LPNR,vNames = dichoVar)
  cohortMean <- round(apply(allData[,hc],2,mean,na.rm=T),1) 
  
  t.u <- rbind(t11.u,t12.u)
  cc <- round(c(cohortMean,cohortPerc),1)
  mm.u <- cbind(cc,t.u)
  cohort <- nrow(allData)
  casesU <- sapply(c.by.histo.LPNR,length)
  N.u <- c(cohort,casesU)
  mm.u <- rbind(N.u,mm.u)
  colnames(mm.u)[1] <- "Cohort"
  rownames(mm.u) <- rNames
  return(mm.u)
  }


