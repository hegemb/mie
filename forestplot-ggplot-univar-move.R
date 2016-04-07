#############################################
## Set paths:
#############################################
setwd("/Users/hegemb/Documents/jobb/TICE/ovarian cancer :: Mie Jareid/scripts")
#setwd("/Users/hegemb/Documents/jobb/TICE/ovarian cancer :: Mie Jareid/ggplot - forestplott fra nett - reproducable table")
wd <- getwd()
setwd("..")
dataPath <- file.path(getwd(), "data") ## folder to save data files
resultPath <- file.path(getwd(),"results2") ## folder to store results
setwd(wd)

#################################
## Source the forest plot function
#################################
source("THE-forest-plot-function.R")

## LOAD GROUP DATA AND P values from csv file
#groupData<-read.csv(file="groupdata.csv",header=T)

## ---------
## Insert my own data: 
groupData <- read.csv2(file=file.path(resultPath,"univariate-cox-for-forest-plot-HR-and-CI.csv"),header=T,sep=";", dec=".")
gr0 <- groupData
Subgroup <- as.character(groupData$group)
Subgroup
## Remove some of the "-uterus" and "-ovarian" in the names:
Subgroup[grep("All.uterus",Subgroup)] <- sapply(Subgroup[grep("All.uterus",Subgroup)],function(x)substr(x,1,nchar(x)-7))
Subgroup[grep("All.ovarial",Subgroup)] <- sapply(Subgroup[grep("All.ovarial",Subgroup)],function(x)substr(x,1,nchar(x)-8))
Subgroup[grep("All.endometrial",Subgroup)] <- sapply(Subgroup[grep("All.endometrial",Subgroup)],function(x)substr(x,1,nchar(x)-7))
Subgroup[grep("All.epithelial",Subgroup)] <- sapply(Subgroup[grep("All.epithelial",Subgroup)],function(x)substr(x,1,nchar(x)-8))
Subgroup <- factor(Subgroup)
groupData <- cbind(groupData[,1],Subgroup,groupData[,2:ncol(groupData)])
head(groupData)
rownames(groupData) <- groupData$Subgroup

## Import p-values from the heterogenity tests:

p_val <- read.csv2(file=file.path(resultPath,"table4-pvalues-for-heterogenity-univariate-Cox.csv"))

#######################################
## Make forest plot for all carcinomas
## Supp figure 1
#######################################
subGroupNames <- c("All.endometrial", "All.epithelial")
Group <- as.factor(rep("All carcinomas",2))
fileName <- "supplementary-figure-1"
p_value <- p_val[p_val$X %in% c("All.endometrial // All.epithelial"),] ## pick out the right rows with p-values
rownames(p_value) <- p_value$X ## Add row names
gr_name <- c("All.endometrial") ##
p_value <- cbind(Group = gr_name, p_value) ## this will be used to color the boxes in the forest plot according to wheter the heterogenity test is significant or not.
makeForestPlot(subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,pp = p_value,hoyde=8)

#######################
## Forest plot for 
## histological subgroups
## Supp figure 2
##########################
hiN <- c("Endometrioid.clearCell","Serous")
subGroupNames <- c(as.vector(sapply(hiN,function(x)paste(x,c("uterus","ovarian"),sep="-"))))
Group <- as.factor(c(rep("Endometroid and clear cell",2),rep(hiN[2],2)))
fileName <- "supplementary-figure-2"
p_value <- p_val[p_val$X %in% c("Endometrioid.clearCell // Endometrioid.clearCell","Serous // Serous"),] ## pick out the right rows with p-values
gr_name <- subGroupNames[c(1,3)] ##
p_value <- cbind(Group = gr_name, p_value[,2:ncol(p_value)]) ## this will be used to color the boxes in the forest plot according to wheter the heterogenity test is significant or not.
p_value
makeForestPlot(subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,pp = p_value,hoyde=12)
#histoName = fileName; subGroupNames = subGroupNames; groupData = groupData;gr=Group;ff=fileName;pp = p_value
#######################
## Forest plot for 
## moved groups
## Supp figure 3
##########################
ss <- c("Endometrioid.clearCell")
subGroupNames <- as.vector(rbind(sapply(ss,function(x)c(paste(c("All.endometrial.plus.ovary","All.epithelial.minus"),x,sep=".")))))
subGroupNames
Group <- as.factor(c(rep("All carcinomas, reclassified",2)))
fileName <- "supplementary-figure-3"
p_value <- p_val[p_val$X %in% c("All.endometrial.plus.ovary.Endometrioid.clearCell // All.epithelial.minus.Endometrioid.clearCell"),] ## pick out the right rows with p-values
gr_name <- subGroupNames[grep("All.endo",subGroupNames)] ##
p_value <- cbind(Group = gr_name, p_value[,2:ncol(p_value)]) ## this will be used to color the boxes in the forest plot according to wheter the heterogenity test is significant or not.
p_value
makeForestPlot(subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,pp = p_value,hoyde=8)
