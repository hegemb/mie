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
groupData <- read.csv2(file=file.path(resultPath,"univariate-cox-for-forest-plot-HR-and-CI2015-09-03.csv"),header=T,sep=";", dec=".")
gr0 <- groupData
# Subgroup <- rep("ovarian",nrow(groupData))
# Subgroup[grep("uterus",groupData$group)] <- "uterus"
# Subgroup[grep("minus",groupData$group)] <- 
Subgroup <- as.character(groupData$group)
Subgroup
## Remove some of the "-uterus" and "-ovarian" in the names:
Subgroup[grep("All.uterus",Subgroup)] <- sapply(Subgroup[grep("All.uterus",Subgroup)],function(x)substr(x,1,nchar(x)-7))
Subgroup[grep("All.ovarial",Subgroup)] <- sapply(Subgroup[grep("All.ovarial",Subgroup)],function(x)substr(x,1,nchar(x)-8))
Subgroup[grep("All.endometrial",Subgroup)] <- sapply(Subgroup[grep("All.endometrial",Subgroup)],function(x)substr(x,1,nchar(x)-7))
Subgroup[grep("All.epithelial",Subgroup)] <- sapply(Subgroup[grep("All.epithelial",Subgroup)],function(x)substr(x,1,nchar(x)-8))
#Subgroup[grep("ovary.uterus",Subgroup)] <- sapply(Subgroup[grep("ovary.uterus",Subgroup)], function(x)substr(x,1,nchar(x)-7))
#Subgroup <- c(rep("uterus",nrow(groupDataU)),rep("ovarial",nrow(groupDataO))) # indicator for subgroup
Subgroup <- factor(Subgroup)
groupData <- cbind(groupData[,1],Subgroup,groupData[,2:ncol(groupData)])
head(groupData)
rownames(groupData) <- groupData$Subgroup

## Import p-values from the heterogenity tests:

p_val <- read.csv2(file=file.path(resultPath,"table-of-pvalues-for-heterogenity-univariate-Cox-with-moved-groups2015-09-10.csv"))

######## -------------------------################
## Small function to call the forest plot. 
## - hName is the name used when writing the pdf-file with plot to file, and to pick out the right subgroups.
## - refName is the name of the reference group, also used to create name of pdf when stored, and to pick out the right subgroups.
## - fileName, again, used to create name of pdf-file. 
# miniCall <- function(hName,refName,fileName){
#   subGroupNames <- c(refName[1],paste0(refName[1],".plus.ovary.",hName),refName[2],paste0(refName[2],".minus.",hName),paste0(hName,"-uterus"),paste0(hName,".moved"),paste0(hName,"-ovarian"))
#   #Group <- as.factor(c(rep("All",4),rep(hName,3)))
#   Group <- as.factor(c(rep("All",3),rep(hName,4))) ## Denne er ikke riktig, men en workaround fordi jeg ikke får til å gruppere dataene riktig. Det skal være c(rep("All",4),rep(hName,3)). 
#   ## Jeg ser at feilen kommer etter at kommandoen "apply(hl_rows,1,function(x)annotation_custom..." i forest plot, og det er antall grupper som ikke blir reversert slik som resten. Jeg velger uansett å la alt være. 
#   p_value <- p_val[grep(hName, p_val$X),]
#   makeForestPlot(histoName = hName, subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,p_value)
# }
# 
# histNavn <- c("Endometrioid","Clear.cell","Endometrioid.clearCell","Serous")
# 
# ########################################
# ## Make forest plot with "all uterus"
# ## and "all ovarial" as reference group:
# ########################################
# lapply(histNavn,function(x)miniCall(hName = x,refName=c("All.uterus","All.ovarial"),fileName = "all"))
# 
# ##########################################
# ## Make forest plot with "all endometrial" 
# ## and "All epithelial" as reference groups
# ###########################################
# lapply(histNavn,function(x)miniCall(hName = x,refName=c("All.endometrial","All.epithelial"),fileName = "endo.epith"))


#######################################
## Make forest plot with non-epithelial 
## and non-endometrial subtypes.
#######################################
subGroupNames <- c("All.uterus", "All.ovarial","All.endometrial", "All.epithelial","non.endometrial", "non.epithelial")
#Group <- as.factor(c(rep("All",subGroupNames))))
Group <- as.factor(c(rep("All",2),rep("All_endometrial/epithelial",2),rep("Non_endometrial/epithelial",2)))
fileName <- "non-epith-non-endo"
p_value <- p_val[p_val$X %in% c("All.uterus // All.ovarial","All.endometrial // All.epithelial"),] ## pick out the right rows with p-values
rownames(p_value) <- p_value$X ## Add row names
p_value <- rbind(p_value[,2:ncol(p_value)],rep(1,ncol(p_value)-1))
gr_name <- c("All.uterus","All.endometrial","non.endometrial") ##
p_value <- cbind(Group = gr_name, p_value) ## this will be used to color the boxes in the forest plot according to wheter the heterogenity test is significant or not.
                 
makeForestPlot(histoName = fileName, subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,pp = p_value)

#######################
## Forest plot for 
## histological subgroups
##########################
hiN <- c("Endometrioid.clearCell","Serous")
subGroupNames <- c("All.endometrial","All.epithelial",as.vector(sapply(hiN,function(x)paste(x,c("uterus","ovarian"),sep="-"))))
Group <- as.factor(c(rep("All",2),rep(hiN[1],2),rep(hiN[2],2)))
fileName <- "subGroups"
p_value <- p_val[p_val$X %in% c("All.endometrial // All.epithelial","Endometrioid.clearCell // Endometrioid.clearCell","Serous // Serous"),] ## pick out the right rows with p-values
gr_name <- subGroupNames[c(1,3,5)] ##
p_value <- cbind(Group = gr_name, p_value[,2:ncol(p_value)]) ## this will be used to color the boxes in the forest plot according to wheter the heterogenity test is significant or not.
p_value
makeForestPlot(histoName = fileName, subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,pp = p_value)

#######################
## Forest plot for 
## moved groups
##########################
ss <- c("Endometrioid.clearCell","Serous")
#subGroupNames <- as.vector(sapply(ss,function(x)c(paste(c("All.endometrial.plus.ovary","All.epithelial.minus"),x,sep="."),paste(x,"moved-uterus",sep="."))))
#subGroupNames <- as.vector(sapply(ss,function(x)c(paste(c("All.endometrial.plus.ovary","All.epithelial.minus"),x,sep="."))))
subGroupNames <- as.vector(rbind(sapply(ss,function(x)c(paste(c("All.endometrial.plus.ovary","All.epithelial.minus"),x,sep="."))),paste(ss,"moved-uterus",sep=".")))

subGroupNames
Group <- as.factor(c(rep(paste(hiN[1],"moved",sep="-"),3),rep(paste(hiN[2],"moved",sep="-"),3)))
fileName <- "moved"
p_value <- p_val[p_val$X %in% c("All.endometrial.plus.ovary.Endometrioid.clearCell // All.epithelial.minus.Endometrioid.clearCell","All.endometrial.plus.ovary.Serous // All.epithelial.minus.Serous"),] ## pick out the right rows with p-values
gr_name <- subGroupNames[grep("All.endo",subGroupNames)] ##
#p_value <- rbind(p_value[1,2:ncol(p_value)],rep(1,ncol(p_value)-1),p_value[2,2:ncol(p_value)],rep(1,ncol(p_value)-1))
#tmp <- subGroupNames[grep("-uterus",subGroupNames)]
#gr_name <- c(gr_name[1],tmp[1],gr_name[2],tmp[2])
p_value <- cbind(Group = gr_name, p_value[,2:ncol(p_value)]) ## this will be used to color the boxes in the forest plot according to wheter the heterogenity test is significant or not.
p_value
makeForestPlot(histoName = fileName, subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,pp = p_value)
 

#######################################
## Make forest plot to compare
## the subgroups of ovarian cancer
## with the main groups of uterus and 
## ovarian cancer - where the subgroup
## is excluded.
########################################
ss <- c("Endometrioid.clearCell","Serous")
subGroupNames <- as.vector(sapply(ss,function(x)c("All.endometrial",paste(c("All.epithelial.minus"),x,sep="."),paste0(x,"-ovarian"))))
p_value <- cbind(Group=subGroupNames[c(1,3,5)],p_val[1:3,2:ncol(p_val)])
p_value 
p_value[,2:ncol(p_value)] <- 1 ## We dont want to use information on p-values in this plot. Put all to one.
p_value
fileName <- "main-groups-compared-to-doubted-groups"
Group <- c(rep(ss[1],3), rep(ss[2],3))
makeForestPlot(histoName = fileName, subGroupNames = subGroupNames, groupData = groupData,gr=Group,ff=fileName,pp = NA)
