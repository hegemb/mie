## This script is used to generate figure 1 in the paper, displaying a forest plot of the mulitvariate analysis with
## covariates TML, BMI and ever smoker (adjusted for breast cancer history of mother).

#############################################
## Set paths:
#############################################
setwd("/Users/hegemb/Documents/jobb/TICE/ovarian cancer :: Mie Jareid/scripts")
wd <- getwd()
setwd("..")
dataPath <- file.path(getwd(), "data") ## folder to save data files
resultPath <- file.path(getwd(),"results2") ## folder to store results
setwd(wd)

## REQUIRED PACKAGES
require(grid)
require(gridExtra)
require(ggplot2)
require(plyr)


## LOAD GROUP DATA AND P values from csv file
#groupData<-read.csv(file="groupdata.csv",header=T)

## ---------
## Insert my own data: 
groupDataU <- read.csv2(file=file.path(resultPath,"table-of-multivariate-Cox-regression-uterus-parity-breastFeed.csv"),header=T)
levels(groupDataU$X)[levels(groupDataU$X)=="All.endometrial"] <- "All carcinomas"
levels(groupDataU$X)[levels(groupDataU$X)=="Endometrioid.clearCell"] <- "Endometrioid and clear cell"
levels(groupDataU$X)[levels(groupDataU$X)=="All.endometrial.plus.ovary.Endometrioid.clearCell"] <- "All carcinomas, reclassified"
groupDataO <- read.csv2(file=file.path(resultPath,"table-of-multivariate-Cox-regression-ovarian-parity-breastFeed.csv"),header=T)
levels(groupDataO$X)[levels(groupDataO$X)=="All.epithelial"] <- "All carcinomas"
levels(groupDataO$X)[levels(groupDataO$X)=="Endometrioid.clearCell"] <- "Endometrioid and clear cell"
levels(groupDataO$X)[levels(groupDataO$X)=="All.epithelial.minus.Endometrioid.clearCell"] <- "All carcinomas, reclassified"
groupData <- rbind(groupDataU, groupDataO) ## merge
groupData
Subgroup <- c(paste0(rep("Uterine",nrow(groupDataU))," (n=",groupDataU$nEvent,")"),paste0(rep("Ovarian",nrow(groupDataO))," (n=",groupDataO$nEvent,")")) # indicator for subgroup
ID <- c(seq(1,nrow(groupData),by=2),seq(2,nrow(groupData),by=2)) ## ID to make the right order of subtypes
groupData <- cbind(ID,groupData[,1],Subgroup,groupData[,2:ncol(groupData)])
colnames(groupData)[2] <- "Group"
groupData

## Import p-values from heterogenity test:
p_hetero <- read.csv2(file.path(resultPath,"table5-pvalues-for-heterogenity-multivariate-Cox.csv"))

## Make a function to generate the forest plots:
skogsplott <- function(file.name, hoyde, groupData, p_hetero){
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
  CI_Data_royk <- makeSubData("royk_groupever")
  CI_Data_TML <- makeSubData("totalMensLifeBfPar")
  
  groupDataAll <- groupData
  groupData <- groupData[,1:4]
  groupData <- cbind(groupData,P_S=rep(0.1,nrow(groupData)), P_G=rep(.88,nrow(groupData)))
  #### -------
  
  ## SYNTHESIZE SOME PLOT DATA - you can load csv instead
  ## EXPECTS 2 columns - integer for 'ID' matching groupdatacsv
  ## AND 'HR' Hazard Rate
  hazardData<-expand.grid(ID=1:nrow(groupData),HR=1:nrow(groupData))
  hazardData$HR<-1.3-runif(nrow(hazardData))*0.7
  hazardData<-rbind(hazardData,ddply(groupData,.(Group),summarize,ID=max(ID)+0.1,HR=NA)[,2:3])
  hazardData<-rbind(hazardData,data.frame(ID=c(0,-1:(-2-blankRows),max(groupData$ID)+1,max(groupData$ID)+2),HR=NA))
  #CI_Data<-ddply(hazardData[!is.na(hazardData$HR),],.(ID),summarize,low=min(HR),high=max(HR),target=mean(HR))
  
  ## P-value for heterogenity:
  ## MAKE COLOR:
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
                       lab=c(""))
  
  LfLabels<-data.frame(x=rep(15,times=1),
                       y=c(0.5,3.7),
                       lab=c(""))
  
  LegendLabels<-data.frame(x=2.8, ##position in y-direction of final plot.
                           y=c(1), ## position in x-diretion of final plot.
                           lab="BMI") ## will be replaced by correct name subsequently.
  
  
  ## BASIC PLOT
  haz<-ggplot(hazardData,aes(factor(ID),HR))+ labs(x=NULL, y=NULL)
  
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
    geom_segment(data=hoyde,aes(x=4, y=1, xend=ee[1],yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
    geom_point(data=CI_Data_BMI,aes(x = rev(factor(ID)), y = target),shape=22,size=5,fill=ff,vjust=0) + ## Alle boksene. 
    geom_errorbar(data=CI_Data_BMI,aes(x=rev(factor(ID)),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
    
    #scale_y_log10() + 
    #coord_flip() +
    coord_flip(ylim=c(.85,1.2))+ 
    geom_text(data=scaledata,aes(3.1,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
    geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
    geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
    geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
    geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
    geom_segment(aes(x = 3.5, y = .9, xend = 3.5, yend = 1.15)) + ## Horizontal rett strek, markerer tallinja.
    theme_bare
  
  subPlot_BMI
  
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
    geom_segment(data=hoyde,aes(x=4, y=1, xend=ee[1],yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
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
    geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
    geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
    geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
    geom_segment(aes(x = 3.5, y = 0.3, xend = 3.5, yend = 1.5)) + ## Horizontal rett strek, markerer tallinja.
    theme_bare
  
  subPlot_royk
  
  #################
  ## Create subplot 
  ## for TML
  #################
  rr <- range(CI_Data_TML[,2:4])
  rr
  scaledata$HR <- round(seq(from=0.9,to=1.05,length.out=4),2)
  LegendLabels$lab <- "TML" ## rename so that the right name of the variable is displayed in the plot.
  ff <- colorPlot$totalMensLifeBfPar
  
  subPlot_TML <- haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    geom_segment(aes(x = 3.2, y = 1, xend = 3.7, yend = 1)) + ## vertikal strek på HR=1 på tallinja.
    geom_segment(data=hoyde,aes(x=4, y=1, xend=ee[1],yend=1),alpha=0.3,linetype=2,size=0.2)+ ## dottet strek ved HR=1 gjennom alle boxplot.
    #geom_hline(aes(yintercept=1),linetype=2, size=0.5, alpha=.5)+ ## dottet strek ved HR=1 gjennom alle boxplot.
    #geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+
    
    geom_point(data=CI_Data_TML,aes(x = rev(factor(ID)), y = target),shape=22,size=5,fill=ff,vjust=0) + ## Alle boksene. 
    geom_errorbar(data=CI_Data_TML,aes(x=rev(factor(ID)),y=target,ymin =low, ymax=high),width=0.5)+ ## Error-bars, mellom 1Q og 3Q.
    
    #scale_y_log10() + 
    #coord_flip() +
    coord_flip(ylim=c(.85,1.07))+ #
    # coord_flip(ylim=c(ylimit[1]-(ylimit[1]*.5), ylimit[2]+(ylimit[1]*.5)))+ ## Bredden på plottene. Legg til/trekk fra 30%. 
    geom_text(data=scaledata,aes(3.1,HR,label=HR), vjust=0.5, size=dataSize) + ## Tallene på linja under plottet.
    geom_text(data=RtLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) + ## "Hazard Ratio", skrevet over plottet.
    geom_text(data=LegendLabels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) + ## Navn på variabel, under plottet.
    geom_point(data=scaledata,aes(3.5,HR),shape=3,size=3) + ## markerer tall på linja med en strek.
    geom_point(aes(2,12),shape=3,alpha=0,vjust=0) + 
    geom_segment(aes(x = 3.5, y = 0.9, xend = 3.5, yend = 1.15)) + ## Horizontal rett strek, markerer tallinja.
    theme_bare
  
  subPlot_TML
  
  ################################
  ## LEFT PANEL WITH NORMAL SCALE
  ################################
  leftPanel<-haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    coord_flip(ylim=c(0,5.2)) +
    geom_point(aes(x=rev(factor(ID)),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=group_p,aes(rev(factor(y)),0.5,label=Group, fontface="bold"),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=groupData,aes(rev(factor(ID)),1,label=Subgroup),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=LfLabels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, hjust=0, size=titleSize) +
    theme_bare
  
  # Display plot:
  grid.arrange(leftPanel,subPlot_TML,subPlot_BMI,subPlot_royk,widths=rep(1,4), ncol=4, nrow=1)
  
  ## print to file:
  pdf(file.path(resultPath,paste0(file.name,"-multivariate-cox-forest-plot.pdf")),width=12)
  grid.arrange(leftPanel,subPlot_TML,subPlot_BMI,subPlot_royk,widths=rep(1,4), ncol=4, nrow=1)
  dev.off()
}



############
## Generate figure 1: 
skogsplott(file.name = "figure1", hoyde = data.frame(ee=14), p_hetero = p_hetero[1:3,], groupData = groupData[c(1:3,5:7),])
## Generate figure 3:
grData <- groupData[c(4,8),]
grData$ID <- rownames(grData) <- c(1,2)
skogsplott(file.name = "figure3", hoyde = data.frame(ee=8), p_hetero = p_hetero[4,], groupData = grData)
