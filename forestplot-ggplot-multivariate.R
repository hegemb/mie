#############################################
## Set paths:
#############################################
setwd("/Users/hegemb/Documents/jobb/TICE/ovarian cancer :: Mie Jareid/scripts")
#setwd("/Users/hegemb/Documents/jobb/TICE/ovarian cancer :: Mie Jareid/ggplot - forestplott fra nett - reproducable table")
wd <- getwd()
setwd("..")
dataPath <- file.path(getwd(), "data") ## folder to save data files
resultPath <- file.path(getwd(),"results") ## folder to store results
setwd(wd)

## REQUIRED PACKAGES
require(grid)
require(gridExtra)
require(ggplot2)
require(plyr)

############################################
### CUSTOMIZE APPEARANCE WITH THESE     ####
############################################
blankRows<-2   # blank rows under boxplot
titleSize<-4
dataSize<-4
boxColor<-"pink"
############################################
############################################
## BASIC THEMES (SO TO PLOT BLANK GRID)
theme_grid <- theme(
  axis.line = element_blank(), 
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(), 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  axis.ticks.length = unit(0.0001, "mm"),
  axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
  legend.position = "none", 
  panel.background = element_rect(fill = "transparent"), 
  panel.border = element_blank(), 
  panel.grid.major = element_line(colour="grey"), 
  panel.grid.minor = element_line(colour="grey"), 
  panel.margin = unit(c(-0.1,-0.1,-0.1,-0.1), "mm"), 
  plot.margin = unit(c(5,0,5,0.01), "mm")
)

theme_bare <- theme_grid +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

## LOAD GROUP DATA AND P values from csv file
#groupData<-read.csv(file="groupdata.csv",header=T)

## ---------
## Insert my own data: 
groupDataU <- read.csv2(file=file.path(resultPath,"table-of-multivariate-Cox-regression-uterus-parity-breastFeed.csv"),header=T)
levels(groupDataU$X)[levels(groupDataU$X)=="All.endometrial"] <- "All"
groupDataO <- read.csv2(file=file.path(resultPath,"table-of-multivariate-Cox-regression-ovarian-parity-breastFeed.csv"),header=T)
levels(groupDataO$X)[levels(groupDataO$X)=="All.epithelial"] <- "All"
groupData <- rbind(groupDataU, groupDataO) ## merge
groupData
Subgroup <- c(rep("uterus",nrow(groupDataU)),rep("ovarial",nrow(groupDataO))) # indicator for subgroup
ID <- c(seq(1,nrow(groupData),by=2),seq(2,nrow(groupData),by=2)) ## ID to make the right order of subtypes
groupData <- cbind(ID,groupData[,1],Subgroup,groupData[,2:ncol(groupData)])
colnames(groupData)[2] <- "Group"
groupData

groupData <- groupData[order(groupData$ID),] ## Let data get the right order 
#groupData <- groupData[nrow(groupData):1,] ## flip data to get the right order in the plot.
rownames(groupData) <- groupData$ID
groupData
names(groupData)[4] <- "NoP"
#groupData$ID <- as.factor(groupData$ID)

## MAke a function that will split the data in "groupData" into data.frames that can be used in the plot:
makeSubData <- function(navn,grData=groupData){
  tt <- paste(navn,c("LCI","UCI","HR"),sep=".")
  uu <- data.frame(ID=grData$ID,grData[tt])
  names(uu) <- c("ID","low","high","target")
  return(uu)
}
## Apply above function: 
CI_Data_BMI <- makeSubData("BMI")
CI_Data_MOR <- makeSubData("MORYes")
CI_Data_royk <- makeSubData("royk_groupever")
CI_Data_TML <- makeSubData("totalMensLifeBfPar")

groupDataAll <- groupData
groupData <- groupData[,1:4]
groupData <- cbind(groupData,P_S=rep(0.1,6), P_G=rep(.88,6))
#### -------

## SYNTHESIZE SOME PLOT DATA - you can load csv instead
## EXPECTS 2 columns - integer for 'ID' matching groupdatacsv
## AND 'HR' Hazard Rate
hazardData<-expand.grid(ID=1:nrow(groupData),HR=1:6)
hazardData$HR<-1.3-runif(nrow(hazardData))*0.7
hazardData<-rbind(hazardData,ddply(groupData,.(Group),summarize,ID=max(ID)+0.1,HR=NA)[,2:3])
hazardData<-rbind(hazardData,data.frame(ID=c(0,-1:(-2-blankRows),max(groupData$ID)+1,max(groupData$ID)+2),HR=NA))
#CI_Data<-ddply(hazardData[!is.na(hazardData$HR),],.(ID),summarize,low=min(HR),high=max(HR),target=mean(HR))

## P-value for heterogenity:
## MAKE COLOR:
## Import data:
p_hetero <- read.csv2(file.path(resultPath,"table-of-pvalues-for-heterogenity-multivariate-Cox-BfPar2015-08-24.csv"))
levels(p_hetero$X)[levels(p_hetero$X)=="All.endometrial"] <- "All" ## change name of level.
colnames(p_hetero)[1] <- "Group"
nrOfGroups <- nrow(groupData)

whichSign <- apply(p_hetero[,-1],2,function(x)seq(1,nrOfGroups,by=2)[which(x < 0.05)])
whichSign <- lapply(whichSign,function(x)sort(c(x,x+1)))
farge <- rep("white",nrOfGroups)
colorPlot <- lapply(whichSign, function(x) {farge[x] <- "black";return(farge)})



## Make the min/max mean labels
hrlabels<-ddply(hazardData[!is.na(hazardData$HR),],.(ID),summarize,lab=paste(round(mean(HR),2)," (",round(min(HR),2),"-",round(max(HR),2),")",sep=""))

## Points to plot on the log scale
scaledata<-data.frame(ID=0,HR=c(0.6,0.8,1.2,1.8))
#scaledata<-data.frame(ID=0,HR=c(0.8,1.2))

## Pull out the Groups & P values
group_p<-ddply(groupData,.(Group),summarize,P=mean(P_G),y=max(ID)+0.1)
levels(group_p$Group)[levels(group_p$Group)=="Endometroid.clearCell"] <- "Endometroid and clear cell"

## identify the rows to be highlighted, and 
## build a function to add the layers
hl_rows<-data.frame(ID=(1:floor(length(unique(hazardData$ID[which(hazardData$ID>0)]))/2))*2,col="lightgrey")
hl_rows$ID<-hl_rows$ID+blankRows+1
hl_rect<-function(col="white",alpha=0.5){
  rectGrob(   x = 0, y = 0, width = 1, height = 1, just = c("left","bottom"), gp=gpar(alpha=alpha, fill=col))
}

## DATA FOR TEXT LABELS
RtLabels<-data.frame(x=rep(15,times=1),
                     y=c(1),
                     lab=c("Hazard Ratio\n(95% CI)"))

LfLabels<-data.frame(x=rep(15,times=2),
                     y=c(0.5,3.7),
                     lab=c("Subgroup","No. of\nPatients"))

LegendLabels<-data.frame(x=2.8, ##position in y-direction of final plot.
                         y=c(1), ## position in x-diretion of final plot.
                         lab="BMI") ## will be replaced by correct name subsequently.


## BASIC PLOT
haz<-ggplot(hazardData,aes(factor(ID),HR))+ labs(x=NULL, y=NULL)

# ## RIGHT PANEL WITH LOG SCALE
# rightPanel<-haz + 
#   apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
#   geom_segment(aes(x = 2.2, y = 1, xend = 2.7, yend = 1)) + ## vertikal strek på HR=1 på tallinja.
#   geom_hline(aes(yintercept=1),linetype=2, size=0.5)+ ## dottet strek ved HR=1 gjennom alle boxplot.
#   #geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+
#   
#   geom_point(data=CI_Data_BMI,aes(x = factor(ID), y = target),shape=22,size=5,fill=boxColor,vjust=0) + ## Alle boksene. 
#   geom_errorbar(data=CI_Data_BMI,aes(x=factor(ID),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
#   
#   scale_y_log10() + coord_flip() +
#   geom_text(data=scaledata,aes(3,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
#   geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
#  # geom_text(data=hrlabels,aes(factor(ID),4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
#  # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
#  # geom_text(data=groupData,aes(factor(ID),6.5,label=P_S),vjust=0.5, hjust=1, size=dataSize) +
#   geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
#   geom_point(data=scaledata,aes(2.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
#   geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
#   geom_segment(aes(x = 2.5, y = 0.5, xend = 2.5, yend = 3)) + ## Horizontal rett strek, markerer tallinja.
# #  geom_segment(aes(x = 2, y = 1, xend = 2, yend = 1.8),arrow=arrow(),linetype=1,size=1) + 
#  # geom_segment(aes(x = 2, y = 1, xend = 2, yend = 0.6),arrow=arrow(),linetype=1,size=1) + 
#   theme_bare
# 
# rightPanel

# sub_Panel<- function(ssData,sname){
#   trunk <- 3 ## Truncate values at this level, to avoid too many long viskers.
#   tmp_ssData <- ssData[,2:4]
#   tmp_ssData[which(tmp_ssData>3,arr.ind=TRUE)] <- trunk
#   ssData[,2:4] <- tmp_ssData
#   
#   ylimit <- range(ssData[,2:4])
#   LegendLabels$lab <- sname ## rename so that the right name of the variable is displayed in the plot.
#   scaledata$HR <- round(seq(from=ylimit[1],to=ylimit[2],length.out=4),2)
#   subPlot <- haz + 
#     apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
#     geom_segment(aes(x = 3.2, y = 1, xend = 3.7, yend = 1)) + ## vertikal strek på HR=1 på tallinja.
#     geom_segment(aes(x=4, y=1, xend=14,yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
#     #geom_hline(aes(yintercept=1),linetype=2, size=0.5, alpha=.5)+ ## dottet strek ved HR=1 gjennom alle boxplot.
#     #geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+
#     
#     geom_point(data=ssData,aes(x = factor(ID), y = target),shape=22,size=5,fill=boxColor,vjust=0) + ## Alle boksene. 
#     geom_errorbar(data=ssData,aes(x=factor(ID),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
#     
#     #scale_y_log10() + 
#     #coord_flip() +
#     #coord_flip(ylim=c(.3,4))+ #
#     coord_flip(ylim=c(ylimit[1]-(ylimit[1]*.5), ylimit[2]+(ylimit[1]*.5)))+ ## Bredden på plottene. Legg til/trekk fra 30%. 
#     geom_text(data=scaledata,aes(3.1,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
#     geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
#     # geom_text(data=hrlabels,aes(factor(ID),4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
#     # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
#     # geom_text(data=groupData,aes(factor(ID),6.5,label=P_S),vjust=0.5, hjust=1, size=dataSize) +
#     geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
#     geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
#     geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
#     geom_segment(aes(x = 3.5, y = ylimit[1], xend = 3.5, yend = ylimit[2])) + ## Horizontal rett strek, markerer tallinja.
#     #  geom_segment(aes(x = 2, y = 1, xend = 2, yend = 1.8),arrow=arrow(),linetype=1,size=1) + 
#     # geom_segment(aes(x = 2, y = 1, xend = 2, yend = 0.6),arrow=arrow(),linetype=1,size=1) + 
#     theme_bare
#   return(subPlot)
# }
# 
# ## Apply the funciton above to make panels of plots, that will be placed next to each other:
# #BMIp <- sub_Panel(CI_Data_BMI,"BMI")
# BMIp <- sub_Panel(CI_Data_BMI,"BMI")
# MORp <- sub_Panel(CI_Data_MOR,"Bc. of mother")
# roykp <- sub_Panel(CI_Data_royk,"Smoking history")
# TMLp <- sub_Panel(CI_Data_TML,"Tot.mens.life")

##########################
### CREATE PANEL FOR BMI: 
###########################
#sub_Panel_BMI<- function(ssData,sname){
   ff <- colorPlot$BMI ## Farge på boksen avgjøres av om p-verdi for heterogenityet er signifikant eller ei.
  scaledata$HR <- c(0.95,1,1.05,1.10)#round(seq(from=0.9,to=1.15,length.out=4),2)
LegendLabels$lab <- "BMI" ## rename so that the right name of the variable is displayed in the plot.
  subPlot_BMI <- haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    geom_segment(aes(x = 3.2, y = 1, xend = 3.7, yend = 1)) + ## vertikal strek på HR=1 på tallinja.
    geom_segment(aes(x=4, y=1, xend=14,yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
    #geom_hline(aes(yintercept=1),linetype=2, size=0.5, alpha=.5)+ ## dottet strek ved HR=1 gjennom alle boxplot.
    #geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+
    
    geom_point(data=CI_Data_BMI,aes(x = rev(factor(ID)), y = target),shape=22,size=5,fill=ff,vjust=0) + ## Alle boksene. 
    geom_errorbar(data=CI_Data_BMI,aes(x=rev(factor(ID)),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
    
    #scale_y_log10() + 
    #coord_flip() +
    coord_flip(ylim=c(.85,1.2))+ #
   # coord_flip(ylim=c(ylimit[1]-(ylimit[1]*.5), ylimit[2]+(ylimit[1]*.5)))+ ## Bredden på plottene. Legg til/trekk fra 30%. 
    geom_text(data=scaledata,aes(3.1,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
    geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
    # geom_text(data=hrlabels,aes(factor(ID),4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=groupData,aes(factor(ID),6.5,label=P_S),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
    geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
    geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
    geom_segment(aes(x = 3.5, y = .9, xend = 3.5, yend = 1.15)) + ## Horizontal rett strek, markerer tallinja.
    #  geom_segment(aes(x = 2, y = 1, xend = 2, yend = 1.8),arrow=arrow(),linetype=1,size=1) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = 0.6),arrow=arrow(),linetype=1,size=1) + 
    theme_bare

subPlot_BMI
#################
## Create subplot 
## for MOR 
#################
range(CI_Data_MOR[,2:4])
scaledata$HR <- c(0.5,2,3,4)#round(seq(from=0.3,to=4,length.out=4),1)
LegendLabelsMOR <- LegendLabels
LegendLabelsMOR$lab <- "Bc. hist. of mother" ## rename so that the right name of the variable is displayed in the plot.
LegendLabelsMOR$y <- 1.5
ff <- colorPlot$MORYes
subPlot_MOR <- haz + 
  apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
  geom_segment(aes(x = 3.2, y = 1, xend = 3.7, yend = 1)) + ## vertikal strek på HR=1 på tallinja.
  geom_segment(aes(x=4, y=1, xend=14,yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
  #geom_hline(aes(yintercept=1),linetype=2, size=0.5, alpha=.5)+ ## dottet strek ved HR=1 gjennom alle boxplot.
  #geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+
  
  geom_point(data=CI_Data_MOR,aes(x = rev(factor(ID)), y = target),shape=22,size=5,fill=ff,vjust=0) + ## Alle boksene. 
  geom_errorbar(data=CI_Data_MOR,aes(x=rev(factor(ID)),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
  
  #scale_y_log10() + 
  #coord_flip() +
  coord_flip(ylim=c(.1,5))+ #
  # coord_flip(ylim=c(ylimit[1]-(ylimit[1]*.5), ylimit[2]+(ylimit[1]*.5)))+ ## Bredden på plottene. Legg til/trekk fra 30%. 
  geom_text(data=scaledata,aes(3.1,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
  geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
  # geom_text(data=hrlabels,aes(factor(ID),4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
  # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
  # geom_text(data=groupData,aes(factor(ID),6.5,label=P_S),vjust=0.5, hjust=1, size=dataSize) +
  geom_text(data=LegendLabelsMOR,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
  geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
  geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
  geom_segment(aes(x = 3.5, y = .3, xend = 3.5, yend = 4)) + ## Horizontal rett strek, markerer tallinja.
  #  geom_segment(aes(x = 2, y = 1, xend = 2, yend = 1.8),arrow=arrow(),linetype=1,size=1) + 
  # geom_segment(aes(x = 2, y = 1, xend = 2, yend = 0.6),arrow=arrow(),linetype=1,size=1) + 
  theme_bare

subPlot_MOR

#################
## Create subplot 
## for royk
#################
rr <- range(CI_Data_royk[,2:4])
rr
scaledata$HR <- round(seq(from=rr[1],to=rr[2],length.out=4),1)
LegendLabels$lab <- "Ever smoker" ## rename so that the right name of the variable is displayed in the plot.
ff <- colorPlot$royk_groupever

subPlot_royk <- haz + 
  apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
  geom_segment(aes(x = 3.2, y = 1, xend = 3.7, yend = 1)) + ## vertikal strek på HR=1 på tallinja.
  geom_segment(aes(x=4, y=1, xend=14,yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
  #geom_hline(aes(yintercept=1),linetype=2, size=0.5, alpha=.5)+ ## dottet strek ved HR=1 gjennom alle boxplot.
  #geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+
  
  geom_point(data=CI_Data_royk,aes(x = rev(factor(ID)), y = target),shape=22,size=5,fill=ff,vjust=0) + ## Alle boksene. 
  geom_errorbar(data=CI_Data_royk,aes(x=rev(factor(ID)),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
  
  #scale_y_log10() + 
  #coord_flip() +
  coord_flip(ylim=c(.2,1.7))+ #
  # coord_flip(ylim=c(ylimit[1]-(ylimit[1]*.5), ylimit[2]+(ylimit[1]*.5)))+ ## Bredden på plottene. Legg til/trekk fra 30%. 
  geom_text(data=scaledata,aes(3.1,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
  geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
  # geom_text(data=hrlabels,aes(factor(ID),4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
  # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
  # geom_text(data=groupData,aes(factor(ID),6.5,label=P_S),vjust=0.5, hjust=1, size=dataSize) +
  geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
  geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
  geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
  geom_segment(aes(x = 3.5, y = 0.3, xend = 3.5, yend = 1.5)) + ## Horizontal rett strek, markerer tallinja.
  #  geom_segment(aes(x = 2, y = 1, xend = 2, yend = 1.8),arrow=arrow(),linetype=1,size=1) + 
  # geom_segment(aes(x = 2, y = 1, xend = 2, yend = 0.6),arrow=arrow(),linetype=1,size=1) + 
  theme_bare

subPlot_royk

#################
## Create subplot 
## for TML
#################
rr <- range(CI_Data_TML[,2:4])
rr
scaledata$HR <- c(0.95,1,1.05,1.10)#round(seq(from=0.95,to=1.15,length.out=4),2)
LegendLabels$lab <- "Tot. mens. lifespan" ## rename so that the right name of the variable is displayed in the plot.
ff <- colorPlot$totalMensLifeBfPar

subPlot_TML <- haz + 
  apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
  geom_segment(aes(x = 3.2, y = 1, xend = 3.7, yend = 1)) + ## vertikal strek på HR=1 på tallinja.
  geom_segment(aes(x=4, y=1, xend=14,yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
  #geom_hline(aes(yintercept=1),linetype=2, size=0.5, alpha=.5)+ ## dottet strek ved HR=1 gjennom alle boxplot.
  #geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+
  
  geom_point(data=CI_Data_TML,aes(x = rev(factor(ID)), y = target),shape=22,size=5,fill=ff,vjust=0) + ## Alle boksene. 
  geom_errorbar(data=CI_Data_TML,aes(x=rev(factor(ID)),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
  
  #scale_y_log10() + 
  #coord_flip() +
  coord_flip(ylim=c(.9,1.2))+ #
  # coord_flip(ylim=c(ylimit[1]-(ylimit[1]*.5), ylimit[2]+(ylimit[1]*.5)))+ ## Bredden på plottene. Legg til/trekk fra 30%. 
  geom_text(data=scaledata,aes(3.1,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
  geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
  # geom_text(data=hrlabels,aes(factor(ID),4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
  # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
  # geom_text(data=groupData,aes(factor(ID),6.5,label=P_S),vjust=0.5, hjust=1, size=dataSize) +
  geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
  geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
  geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
  geom_segment(aes(x = 3.5, y = 0.95, xend = 3.5, yend = 1.15)) + ## Horizontal rett strek, markerer tallinja.
  #  geom_segment(aes(x = 2, y = 1, xend = 2, yend = 1.8),arrow=arrow(),linetype=1,size=1) + 
  # geom_segment(aes(x = 2, y = 1, xend = 2, yend = 0.6),arrow=arrow(),linetype=1,size=1) + 
  theme_bare

subPlot_TML

# boxplotForest <- function(colNavn, labNavn, haz,ylimit){
#   hData <- groupDataAll[,c("ID",paste(colNavn,c("HR","LCI","UCI"),sep="."))]
#   colnames(hData) <- c("ID","HR","LCI","UCI")
#   CI_Data<-data.frame(ID=as.numeric(hData$ID),low=hData$LCI,high=hData$UCI,target=hData$HR)
#   # hazardData <- hData
#   # rownames(hazardData) <- hazardData$ID
#   LegendLabels<-data.frame(x=0.7,y=1,lab=labNavn)
#   #scaledata <- data.frame(ID=rep(0,5),HR=quantile(as.matrix(CI_Data[,2:4])))
#   ff <- colorPlot[[colNavn]]  
#   
#   boxPanel<-haz + 
#     apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
#     coord_flip(ylim=c(ylimit))+ # ylim her justerer bredden på plottene.
#     #geom_segment(aes(x = 2, y = 1, xend = 1.5, yend = 1)) + 
#     #geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
#     ## Boxplots
#     geom_hline(aes(yintercept=1),linetype=2, size=0.5) + 
#     geom_point(data=CI_Data,aes(x = rev(factor(ID)), y = target),shape=21,size=3,fill=ff,vjust=0) + ## her justeres størrelse på HR punkt og form.
#     geom_errorbar(data=CI_Data,aes(x=rev(factor(ID)),y=target,ymin =low, ymax=high),width=0.5) + 
#     #geom_text(data=RtLabels,aes(x[1],y[1],label=lab[1], fontface="bold"), vjust=1, size=titleSize) +
#     geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.3, vjust=1, size=titleSize) + 
#     #geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"), size=titleSize) +   
#     #geom_text(data=scaledata,aes(0,HR,label=HR), vjust=0.5, size=dataSize) +
#     #geom_point(aes(2,12),shape=3,alpha=0,vjust=0) +
#     #     geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 13)) + 
#     #     geom_segment(aes(x = 2, y = 1, xend = 2, yend = 1.8),arrow=arrow(),linetype=1,size=1) + 
#     #     geom_segment(aes(x = 2, y = 1, xend = 2, yend = 0.2),arrow=arrow(),linetype=1,size=1) + 
#     #     #scale_y_continous("BMI", limits=c(-2,2)) +
#     theme_bare
#   return(boxPanel)
# }
# 
# 
# TMLpanel <- boxplotForest(colNavn = "totalMensLifeBfPar", labNavn= "Tot.Mens.Life",haz,ylimit=c(0.9,1.4))


## LEFT PANEL WITH NORMAL SCALE
leftPanel<-haz + 
  apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
  coord_flip(ylim=c(0,5.2)) +
  geom_point(aes(x=rev(factor(ID)),y=1),shape=3,alpha=0,vjust=0) + 
  geom_text(data=group_p,aes(rev(factor(y)),0.5,label=Group, fontface="bold"),vjust=0.5, hjust=0, size=dataSize) +
  geom_text(data=groupData,aes(rev(factor(ID)),1,label=Subgroup),vjust=0.5, hjust=0, size=dataSize) +
  geom_text(data=groupData,aes(rev(factor(ID)),5,label=NoP),vjust=0.5, hjust=1, size=dataSize) +
  geom_text(data=LfLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, hjust=0, size=titleSize) +
  #geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + ## linje under plottet - "tallinja".
  theme_bare

#grid.arrange(leftPanel,BMIp, widths=rep(1.5,2), ncol=2, nrow=1)
## PLOT THEM BOTH IN A GRID SO THEY MATCH UP
# pdf(file.path(resultPath,paste0("tubal-multivariate-cox-forest-plot",Sys.Date(),".pdf")),width=12)
# grid.arrange(leftPanel,BMIp,MORp,roykp,TMLp, widths=rep(1,5), ncol=5, nrow=1)
# dev.off()

## Display plot:
grid.arrange(leftPanel,subPlot_TML,subPlot_BMI,subPlot_royk,subPlot_MOR, widths=rep(1,5), ncol=5, nrow=1)

## print to file:
pdf(file.path(resultPath,paste0("tubal-multivariate-cox-forest-plot",Sys.Date(),".pdf")),width=12)
grid.arrange(leftPanel,subPlot_TML,subPlot_BMI,subPlot_royk,subPlot_MOR, widths=rep(1,5), ncol=5, nrow=1)
dev.off()