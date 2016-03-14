## This script holds a function that again will call 
## other functions in order to generate
## output needed for Cox ph relative risk calculations, 
## checking the Cox ph assumptions and generate the 
## output as readable plots/tables.
## The functions that we import and the workflow is the 
## same as that descripbed in "main-ovarian.R", but made "automatic".
## The reason is that we need to be able to easly generate output 
## for different subsets of the data.

## The input of this function must be 
## - variable.list: A matrix (nxp) of p covariates for n individuals. 
## - startTid: Vector of length n with inclusion in study times. 
## - stoppTid: Vector of length n with end of follow up.
## - event: Vector of length n with event indicators (0 - censoring, 1 - event). 
## - navn: Name used to generate files.
## The function returns a coxph object.

do_all_cox <- function(data.sub, navn){
  
  #############################################
  ## Import libraries:
  #############################################
  require(survival)
  require(knitr)
  #install.packages("pander")
  #require(pander)
  
  ###########################
  ## Import functions:
  ###########################
  source("do_calc_relative_risk.R") ## Relative risk. 
  #source("do_diagnose_propHazard.R") ## Check Cox PH assumptions. 
  source("do_make_table_HR.R") ## Produce output.
  
  ###############################
  ## Calculate relative risk,  ##
  ## check cox assumptions,    ##
  ## and make table of output: ##
  ###############################
  
  ## CALCULATE RELATIVE RISK: ##
  ## Extract variables to use in the Cox regression:
  lpnr <- data.sub[,1]
  startAld <- data.sub[,2]
  stoppAld <- data.sub[,3]
  event <- data.sub[,4]
  #print(paste("\n\n\n",name,"\n Number of events",sum(event)))
  variable.list <- data.sub[,5:ncol(data.sub)]
  
  ## Apply Cox regression with left truncated and right censored data:
  coxfit <- try(do_calc_relative_risk(variable.list, startAld, stoppAld, event)) ## Returns a list of coxph objects
  
  ## CHECK COX ASSUMPTIONS ## 
  ## Make plots:
  #   print("Check Cox assumptions and save plots to file")
  #   vvName <- names(coxfit) ## extract variable names (equal for all cancer types).
  #   pdf(file.path(resultPath,paste0("check-cox-assumptions-",navn,".pdf")))
  #   for(j in 1:length(vvName)){
  #     #print(paste(vvName[j], navn, sep="-"))
  #     do_diagnose_propHazard(x=coxfit[[j]],cancerName=navn,varName=vvName[j])
  #   }
  #   dev.off()
  #   
  ## MAKE TABLE OF OUTPUT ##
  outTab <- do_make_table_HR(coxfit, variable.list, event, navn, lpnr)
  
  return(coxfit)
}
