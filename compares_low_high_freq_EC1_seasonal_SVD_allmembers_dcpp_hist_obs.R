# This script compares the Expansion Coefficients of SVD1 of SST/SLP for DCPP,
# Historical & Observations
#
# M. Alvarez (2021)
#--------------------------------------------------------------------------------

# Call libraries to be used
library("ncdf4")
library("metR")
library('pracma')
library('lubridate')
library('reshape2')
library('data.table')
library("RColorBrewer")
library("maps")
library("ggplot2")
library("gridExtra")
library("grid")
library("maptools")
library("directlabels")
library("zoo")
library("ggnewscale")

rm(list=ls())
graphics.off()

setwd("/home/maralv/")

#---------------------------------------------------------------------------------------
#  functions 
#---------------------------------------------------------------------------------------
scaleFUN <- function(x) sprintf("%.1f", x)

#---------------------------------------------------------------------------------------
#  Man Program
#---------------------------------------------------------------------------------------
#################### Settings ########################
# # Remove trend? TRUE or FALSE
remove.trend=TRUE
# Select season
sel.season="JJA"
# Lead selection inside a for loop
######################################################

# Data from major volcanic eruptions since 1961
volcanic = data.table(date = c(as.Date("1963-03-17"),as.Date("1982-03-28"),as.Date("1991-06-15")),name=c("Agung","El Chichón","Pinatubo"),plot.value=-c(1.85,1.85,1.85))
# Data from EL Ninio (ONI index). Threshold: +/- 0.5
enso = as.data.table(read.table("/home/maralv/Documents/ONI.Nino.Index.csv",header=TRUE))

if(sel.season == "DJF"){
  # Adjust ENSO values
  enso = enso[,.(Year,DJF)]
  enso = enso[,date := as.Date(paste0(as.character(Year),"-01-16"))]
  enso$plot.value = 2.75
  setnames(enso,"DJF","oni")
  
}else if(sel.season == "JJA"){
  # Adjust ENSO values
  enso = enso[,.(Year,JJA)]
  enso = enso[,date := as.Date(paste0(as.character(Year),"-07-16"))]
  enso$plot.value = 2.75
  setnames(enso,"JJA","oni")
}

# Create data table to save correlation coefficients
# SST
corr.dcpp.sst.lf = data.table(member=1:16)
corr.dcpp.sst.lf[member==16]$member=90
corr.dcpp.sst.lf = corr.dcpp.sst.lf[, .(lead = 1:10), by = c(colnames(corr.dcpp.sst.lf))]
corr.dcpp.sst.lf$r = NA_real_


corr.dcpp.sst.hf = data.table(member=1:16)
corr.dcpp.sst.hf[member==16]$member=90
corr.dcpp.sst.hf = corr.dcpp.sst.hf[, .(lead = 1:10), by = c(colnames(corr.dcpp.sst.hf))]
corr.dcpp.sst.hf$r = NA_real_

# SLP
corr.dcpp.slp.lf = data.table(member=1:16)
corr.dcpp.slp.lf[member==16]$member=90
corr.dcpp.slp.lf = corr.dcpp.slp.lf[, .(lead = 1:10), by = c(colnames(corr.dcpp.slp.lf))]
corr.dcpp.slp.lf$r = NA_real_


corr.dcpp.slp.hf = data.table(member=1:16)
corr.dcpp.slp.hf[member==16]$member=90
corr.dcpp.slp.hf = corr.dcpp.slp.hf[, .(lead = 1:10), by = c(colnames(corr.dcpp.slp.hf))]
corr.dcpp.slp.hf$r = NA_real_

##############
# # Start for on leads
for(l in 1:10){
sel.lead=l

# sel.lead=1

# Load data

# Observations
load(paste0("/home/maralv/data/ECs_SVD1_",sel.season,"_OBS_weighted_notrend.rda"))
sst.OBS=sst
psl.OBS=slp
exp.coef.norm.OBS=exp.coef.norm
rm(sst,slp,exp.coef.norm)

# Historical
load(paste0("/home/maralv/data/ECs_SVD1_",sel.season,"_historical_weighted_notrend_allmembers.rda"))
sst.ECE.allin1.hist=sst.ECE.allin1
psl.ECE.allin1.hist=psl.ECE.allin1
exp.coef.separated.hist=exp.coef.separated
rm(sst.ECE.allin1,psl.ECE.allin1,exp.coef.separated)

# Decadal Predictions
load(paste0("/home/maralv/data/ECs_SVD1_",sel.season,"_dcpp_weighted_notrend_allmembers_lead_year_",sel.lead,".rda"))
sst.ECE.allin1.dcpp=sst.ECE.allin1
psl.ECE.allin1.dcpp=psl.ECE.allin1
exp.coef.separated.dcpp=exp.coef.separated
rm(sst.ECE.allin1,psl.ECE.allin1,exp.coef.separated)

###### Making targetdates match and falling in the middle of the season
if(sel.season=="DJF"){
  # DCPP: Targetdate year is ok by construction, mmdd set to 01-16
  exp.coef.separated.dcpp$targetdate = as.Date(as.character(paste0(year(exp.coef.separated.dcpp$targetdate),"-01-16")))
  
  # Historical: date set as 1949-12-16 to represent DJF 49/50 -> a year must be added
  exp.coef.separated.hist=exp.coef.separated.hist[date>=as.Date("1949-12-01")]                                  # Eliminates 1949 data because Dec. 1948 was not used
  exp.coef.separated.hist$date = exp.coef.separated.hist$date %m+% years(1)                                     # Add 1 year
  exp.coef.separated.hist$date = as.Date(as.character(paste0(year(exp.coef.separated.hist$date),"-01-16")))     # Moves all dates to 16/01  
  
  # Obs: Date taken from December month (mmdd set to 12-01) -> a year must be added
  exp.coef.norm.OBS=exp.coef.norm.OBS[date>=as.Date("1949-12-01")]                                  # Eliminates 1949 data because Dec. 1948 was not used
  exp.coef.norm.OBS$date = exp.coef.norm.OBS$date %m+% years(1)                                     # Add 1 year
  exp.coef.norm.OBS$date = as.Date(as.character(paste0(year(exp.coef.norm.OBS$date),"-01-16")))     # Moves all dates to 01/01

}else if(sel.season=="JJA"){
  # No problem with dates, mmdd should be aligned.
  
  # DCPP: Targetdate year is ok by construction, mmdd set to 06-16
  exp.coef.separated.dcpp$targetdate = as.Date(as.character(paste0(year(exp.coef.separated.dcpp$targetdate),"-07-16")))
  
  # Historical: date set as 1949-06-16 to represent JJA 49 -> a year must be added
  exp.coef.separated.hist$date = as.Date(as.character(paste0(year(exp.coef.separated.hist$date),"-07-16")))     # Moves all dates to 06/01  
  
  # Obs: Date set to 06-01 -> a year must be added
  exp.coef.norm.OBS$date = as.Date(as.character(paste0(year(exp.coef.norm.OBS$date),"-07-16")))     # Moves all dates to 06/01, just in case
  
}


# Reduce period for comparison (1961-2014)
exp.coef.norm.OBS = exp.coef.norm.OBS[date>=as.Date("1961-01-01") & date<as.Date("2015-01-01")]
exp.coef.separated.dcpp = exp.coef.separated.dcpp[targetdate>=as.Date("1961-01-01") & targetdate<as.Date("2015-01-01")]
exp.coef.separated.hist = exp.coef.separated.hist[date>=as.Date("1961-01-01") & date<as.Date("2015-01-01")]
enso = enso[date>=as.Date("1961-01-01") & date<as.Date("2015-01-01"),]

# Merge OBS & DCPP EC1 time series
setnames(exp.coef.separated.dcpp,"targetdate","date")
exp.coef.norm.OBS = exp.coef.norm.OBS[, .(member = 1:15), by = c(colnames(exp.coef.norm.OBS))]
ec1.obs.dcpp=merge(exp.coef.norm.OBS[,.(date,ec.sst.1,ec.slp.1,member)],exp.coef.separated.dcpp[,.(member,date,ec1.sst.norm,ec1.slp.norm,sst.med,slp.med)],by=c("date","member"),all=TRUE)
setnames(ec1.obs.dcpp,c("ec1.sst.norm","ec1.slp.norm","sst.med","slp.med"),c("ec1.sst.norm.dcpp","ec1.slp.norm.dcpp","sst.med.dcpp","slp.med.dcpp"))

# Compute low frequency EC1s
# OBS
ec1.obs.dcpp=ec1.obs.dcpp[,ec.sst.1.lf := frollmean(ec.sst.1,5,fill=NA,align="center"),by="member"]
ec1.obs.dcpp$ec.sst.1.hf = ec1.obs.dcpp$ec.sst.1 - ec1.obs.dcpp$ec.sst.1.lf
ec1.obs.dcpp=ec1.obs.dcpp[,ec.slp.1.lf := frollmean(ec.slp.1,5,fill=NA,align="center"),by="member"]
ec1.obs.dcpp$ec.slp.1.hf = ec1.obs.dcpp$ec.slp.1 - ec1.obs.dcpp$ec.slp.1.lf
# DCPP
ec1.obs.dcpp=ec1.obs.dcpp[,ec1.sst.norm.dcpp.lf := frollmean(ec1.sst.norm.dcpp,5,fill=NA,align="center"),by="member"]
ec1.obs.dcpp$ec1.sst.norm.dcpp.hf = ec1.obs.dcpp$ec1.sst.norm.dcpp - ec1.obs.dcpp$ec1.sst.norm.dcpp.lf
ec1.obs.dcpp=ec1.obs.dcpp[,ec1.slp.norm.dcpp.lf := frollmean(ec1.slp.norm.dcpp,5,fill=NA,align="center"),by="member"]
ec1.obs.dcpp$ec1.slp.norm.dcpp.hf = ec1.obs.dcpp$ec1.slp.norm.dcpp - ec1.obs.dcpp$ec1.slp.norm.dcpp.lf

# Compute max,p80,median,p20,min by member
ec1.obs.dcpp[,sst.min.lf := min(ec1.sst.norm.dcpp.lf),by=c("date")]
ec1.obs.dcpp[,sst.p20.lf := quantile(ec1.sst.norm.dcpp.lf,0.13333,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,sst.med.lf := mean(ec1.sst.norm.dcpp.lf),by=c("date")]
ec1.obs.dcpp[,sst.p80.lf := quantile(ec1.sst.norm.dcpp.lf,0.86667,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,sst.max.lf := max(ec1.sst.norm.dcpp.lf),by=c("date")]

# Compute max,p80,median,p20,min by member
ec1.obs.dcpp[,slp.min.lf := min(ec1.slp.norm.dcpp.lf),by=c("date")]
ec1.obs.dcpp[,slp.p20.lf := quantile(ec1.slp.norm.dcpp.lf,0.13333,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,slp.med.lf := mean(ec1.slp.norm.dcpp.lf),by=c("date")]
ec1.obs.dcpp[,slp.p80.lf := quantile(ec1.slp.norm.dcpp.lf,0.86667,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,slp.max.lf := max(ec1.slp.norm.dcpp.lf),by=c("date")]


# Compute max,p80,median,p20,min by member
ec1.obs.dcpp[,sst.min.hf := min(ec1.sst.norm.dcpp.hf),by=c("date")]
ec1.obs.dcpp[,sst.p20.hf := quantile(ec1.sst.norm.dcpp.hf,0.13333,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,sst.med.hf := mean(ec1.sst.norm.dcpp.hf),by=c("date")]
ec1.obs.dcpp[,sst.p80.hf := quantile(ec1.sst.norm.dcpp.hf,0.86667,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,sst.max.hf := max(ec1.sst.norm.dcpp.hf),by=c("date")]

# Compute max,p80,median,p20,min by member
ec1.obs.dcpp[,slp.min.hf := min(ec1.slp.norm.dcpp.hf),by=c("date")]
ec1.obs.dcpp[,slp.p20.hf := quantile(ec1.slp.norm.dcpp.hf,0.13333,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,slp.med.hf := mean(ec1.slp.norm.dcpp.hf),by=c("date")]
ec1.obs.dcpp[,slp.p80.hf := quantile(ec1.slp.norm.dcpp.hf,0.86667,na.rm=TRUE),by=c("date")]
ec1.obs.dcpp[,slp.max.hf := max(ec1.slp.norm.dcpp.hf),by=c("date")]


##########
# Plot
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")
breaks.dates=exp.coef.norm.OBS[member==1]$date
limits.EC=2

# # Version all
# g1 <- ggplot() +
#   theme_bw()+
#   geom_line(data=ec1.obs.dcpp,aes(date, ec1.sst.norm.dcpp.lf,group=member,color=as.factor(member)),alpha=0.8,size=0.5)+
#   scale_color_manual(values=pal)+
#   geom_line(data=ec1.obs.dcpp[member==1,],aes(date, ec.sst.1.lf),col="black",alpha=0.8,size=1)+
#   theme(legend.position = "none")+
#   scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
#   theme(text = element_text(size=14))+
#   scale_y_continuous(breaks = seq(-limits.EC,limits.EC,0.5),limits=c(-limits.EC,limits.EC),expand = c(0, 0.))+
#   theme(axis.text.x = element_text(angle = 45, hjust=1))+
#   xlab("date")+ylab("EC1 (normalized units)")+
#   ggtitle(paste0("SST EC1-Low Frequency filtered. DCPP lead ",sel.lead))+
#   theme(axis.title = element_text(size=11),title = element_text(size=10))

# Version Ribbon
 g2 <-  ggplot() +
  theme_bw()+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=sst.min.lf , ymax=sst.max.lf ),fill="#de77ae",alpha=0.2)+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=sst.p20.lf , ymax=sst.p80.lf ),fill="#de77ae",alpha=0.4)+
  geom_line(data=ec1.obs.dcpp[member==1,],aes(date, ec.sst.1.lf),col="black",alpha=0.8,size=1)+
    
  # Volcanos
  # geom_point(data=volcanic,aes((date %m+% months(3)),(plot.value+0.3)),shape=19,size=2,color="#bababa")+  
  # geom_point(data=volcanic,aes(date,(plot.value+0.25)),shape=19,size=1.5,color="#e0e0e0")+
  geom_point(data=volcanic,aes(date,plot.value),shape=24,size=4,fill="#8c510a")+
  geom_point(data=volcanic,aes(date,(plot.value+0.1)),shape=17,size=1.5,color="#f46d43")+
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-limits.EC,limits.EC,0.5),limits=c(-limits.EC,limits.EC),expand = c(0, 0.),labels=scaleFUN)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.minor.y=element_blank())+
  xlab("date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1-Low Frequency filtered. DCPP lead ",sel.lead, ", SST"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))


# SLP
# g3 <- ggplot() +
#   theme_bw()+
#   geom_line(data=ec1.obs.dcpp,aes(date, ec1.slp.norm.dcpp.lf,group=member,color=as.factor(member)),alpha=0.8,size=0.5)+
#   scale_color_manual(values=pal)+
#   geom_line(data=ec1.obs.dcpp[member==1,],aes(date, ec.slp.1.lf),col="black",alpha=0.8,size=1)+
#   theme(legend.position = "none")+
#   scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
#   theme(text = element_text(size=14))+
#   scale_y_continuous(breaks = seq(-limits.EC,limits.EC,0.5),limits=c(-limits.EC,limits.EC),expand = c(0, 0.))+
#   theme(axis.text.x = element_text(angle = 45, hjust=1))+
#   xlab("date")+ylab("EC1 (normalized units)")+
#   ggtitle(paste0("SLP EC1-Low Frequency filtered. DCPP lead ",sel.lead))+
#   theme(axis.title = element_text(size=11),title = element_text(size=10))

# Version Ribbon
g4 <-  ggplot() +
  theme_bw()+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=slp.min.lf , ymax=slp.max.lf ),fill="#7fbc41",alpha=0.2)+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=slp.p20.lf , ymax=slp.p80.lf ),fill="#7fbc41",alpha=0.4)+
  geom_line(data=ec1.obs.dcpp[member==1,],aes(date, ec.slp.1.lf),col="black",alpha=0.8,size=1)+
    
  # Volcanos
  # geom_point(data=volcanic,aes((date %m+% months(3)),(plot.value+0.3)),shape=19,size=2,color="#bababa")+  
  # geom_point(data=volcanic,aes(date,(plot.value+0.25)),shape=19,size=1.5,color="#e0e0e0")+
  geom_point(data=volcanic,aes(date,plot.value),shape=24,size=4,fill="#8c510a")+
  geom_point(data=volcanic,aes(date,(plot.value+0.1)),shape=17,size=1.5,color="#f46d43")+
    
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-limits.EC,limits.EC,0.5),limits=c(-limits.EC,limits.EC),expand = c(0, 0.),labels=scaleFUN)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.minor.y=element_blank())+
  xlab("date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1-Low Frequency filtered. DCPP lead ",sel.lead, ", SLP"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

#-------------------------------------------------------------
# High frequency
# Version Ribbon
limits.EC=3
g5 <-  ggplot() +
  theme_bw()+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=sst.min.hf , ymax=sst.max.hf ),fill="#de77ae",alpha=0.2)+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=sst.p20.hf , ymax=sst.p80.hf ),fill="#de77ae",alpha=0.4)+
  geom_line(data=ec1.obs.dcpp[member==1,],aes(date, ec.sst.1.hf),col="black",alpha=0.8,size=1)+
    
  # ENSO
  geom_point(data=enso[oni>=0.5,],aes(date,plot.value),shape=24,size=2,fill="#d6604d")+
  geom_point(data=enso[oni<=-0.5,],aes(date,plot.value),shape=25,size=2,fill="#4393c3")+
    
  # Volcanos
  # geom_point(data=volcanic,aes((date %m+% months(3)),(plot.value-0.9+0.4)),shape=19,size=2,color="#bababa")+  
  # geom_point(data=volcanic,aes(date,(plot.value-0.9+0.35)),shape=19,size=1.5,color="#e0e0e0")+
  geom_point(data=volcanic,aes(date,(plot.value-0.9)),shape=24,size=4,fill="#8c510a")+
  geom_point(data=volcanic,aes(date,(plot.value-0.9+0.15)),shape=17,size=1.5,color="#f46d43")+ 
    
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-limits.EC,limits.EC,1),limits=c(-limits.EC,limits.EC),expand = c(0, 0.),labels=scaleFUN)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.minor.y=element_blank())+
  xlab("date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1-High Frequency filtered. DCPP lead ",sel.lead, ", SST"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))


# SLP
# Version Ribbon
g6 <- ggplot() +
  theme_bw()+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=slp.min.hf , ymax=slp.max.hf ),fill="#7fbc41",alpha=0.2)+
  geom_ribbon(data=ec1.obs.dcpp[member==1,], aes(date, ymin=slp.p20.hf , ymax=slp.p80.hf ),fill="#7fbc41",alpha=0.4)+
  
  # ENSO
  geom_point(data=enso[oni>=0.5,],aes(date,plot.value),shape=24,size=2,fill="#d6604d")+
  geom_point(data=enso[oni<=-0.5,],aes(date,plot.value),shape=25,size=2,fill="#4393c3")+
  
  # Volcanos
  # geom_point(data=volcanic,aes((date %m+% months(3)),(plot.value-0.9+0.4)),shape=19,size=2,color="#bababa")+  
  # geom_point(data=volcanic,aes(date,(plot.value-0.9+0.35)),shape=19,size=1.5,color="#e0e0e0")+
  geom_point(data=volcanic,aes(date,(plot.value-0.9)),shape=24,size=4,fill="#8c510a")+
  geom_point(data=volcanic,aes(date,(plot.value-0.9+0.15)),shape=17,size=1.5,color="#f46d43")+ 
  
  geom_line(data=ec1.obs.dcpp[member==1,],aes(date, ec.slp.1.hf),col="black",alpha=0.8,size=1)+
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-limits.EC,limits.EC,1),limits=c(-limits.EC,limits.EC),expand = c(0, 0.),labels=scaleFUN)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.minor.y=element_blank())+
  xlab("date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1-High Frequency filtered. DCPP lead ",sel.lead, ", SLP"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))


##########################################

fig <- grid.arrange(g2,g4,g5,g6, layout_matrix=rbind(cbind(1,2),cbind(3,4)),top = textGrob(paste0(sel.season,", SAOD (SVD1 SST-SLP): Obs & EC-Earth3"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_EC1_DCPP_HIST_OBS_low_high_freq_lead_",sel.lead,".png"),plot=fig,width = 10, height = 6)

##########################################

# Compute correlations for variable and frequency

# Low frequency
# Reduce to selection for sst
aux=ec1.obs.dcpp[,.(date,member,ec.sst.1.lf,ec1.sst.norm.dcpp.lf)]
aux= aux[, r := cor(ec.sst.1.lf,ec1.sst.norm.dcpp.lf,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
cormmb=aux$r

aux=unique(ec1.obs.dcpp[,.(date,ec.sst.1.lf,sst.med.lf)])
# Save correlation coefficients
corr.dcpp.sst.lf[lead==sel.lead]$r=c(cormmb,cor(aux$ec.sst.1.lf,aux$sst.med.lf,use="pairwise.complete.obs"))

# High frequency
# Reduce to selection for sst
aux=ec1.obs.dcpp[,.(date,member,ec.sst.1.hf,ec1.sst.norm.dcpp.hf)]
aux= aux[, r := cor(ec.sst.1.hf,ec1.sst.norm.dcpp.hf,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
cormmb=aux$r

aux=unique(ec1.obs.dcpp[,.(date,ec.sst.1.hf,sst.med.hf)])
# Save correlation coefficients
corr.dcpp.sst.hf[lead==sel.lead]$r=c(cormmb,cor(aux$ec.sst.1.hf,aux$sst.med.hf,use="pairwise.complete.obs"))
rm(aux)


# Low frequency
# Reduce to selection for slp
aux=ec1.obs.dcpp[,.(date,member,ec.slp.1.lf,ec1.slp.norm.dcpp.lf)]
aux= aux[, r := cor(ec.slp.1.lf,ec1.slp.norm.dcpp.lf,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
cormmb=aux$r

aux=unique(ec1.obs.dcpp[,.(date,ec.slp.1.lf,slp.med.lf)])
# Save correlation coefficients
corr.dcpp.slp.lf[lead==sel.lead]$r=c(cormmb,cor(aux$ec.slp.1.lf,aux$slp.med.lf,use="pairwise.complete.obs"))

# High frequency
# Reduce to selection for slp
aux=ec1.obs.dcpp[,.(date,member,ec.slp.1.hf,ec1.slp.norm.dcpp.hf)]
aux= aux[, r := cor(ec.slp.1.hf,ec1.slp.norm.dcpp.hf,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
cormmb=aux$r

aux=unique(ec1.obs.dcpp[,.(date,ec.slp.1.hf,slp.med.hf)])
# Save correlation coefficients
corr.dcpp.slp.hf[lead==sel.lead]$r=c(cormmb,cor(aux$ec.slp.1.hf,aux$slp.med.hf,use="pairwise.complete.obs"))
rm(aux)

} # End lead

##################
# Plotting
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")
breaks=c(-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)

h1 <- ggplot()+
  theme_bw()+
  geom_point(data=corr.dcpp.sst.lf[member<=16,],aes(x=lead,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_line(data=corr.dcpp.sst.lf[member==90,],aes(x=lead,y=r))+
  scale_x_continuous(breaks=seq(1,10,1),expand = c(0.02, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = breaks,limits=c(-0.7,0.7),expand = c(0, 0.))+
  xlab("Lead")+ylab("Correlation")+
  ggtitle(paste0("Correlation between OBS & DCPP EC1-SST: Low frequency"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))

h2 <- ggplot()+
  theme_bw()+
  geom_point(data=corr.dcpp.slp.lf[member<=16,],aes(x=lead,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_line(data=corr.dcpp.slp.lf[member==90,],aes(x=lead,y=r))+
  scale_x_continuous(breaks=seq(1,10,1),expand = c(0.02, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = breaks,limits=c(-0.7,0.7),expand = c(0, 0.))+
  xlab("Lead")+ylab("Correlation")+
  ggtitle(paste0("Correlation between OBS & DCPP EC1-SLP: Low frequency"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))

h3 <- ggplot()+
  theme_bw()+
  geom_point(data=corr.dcpp.sst.hf[member<=16,],aes(x=lead,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_line(data=corr.dcpp.sst.hf[member==90,],aes(x=lead,y=r))+
  scale_x_continuous(breaks=seq(1,10,1),expand = c(0.02, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = breaks,limits=c(-0.7,0.7),expand = c(0, 0.))+
  xlab("Lead")+ylab("Correlation")+
  ggtitle(paste0("Correlation between OBS & DCPP EC1-SST: High frequency"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))

h4 <- ggplot()+
  theme_bw()+
  geom_point(data=corr.dcpp.slp.hf[member<=16,],aes(x=lead,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_line(data=corr.dcpp.slp.hf[member==90,],aes(x=lead,y=r))+
  scale_x_continuous(breaks=seq(1,10,1),expand = c(0.02, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = breaks,limits=c(-0.7,0.7),expand = c(0, 0.))+
  xlab("Lead")+ylab("Correlation")+
  ggtitle(paste0("Correlation between OBS & DCPP EC1-SLP: High frequency"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))

fig <- grid.arrange(h1,h2,h3,h4, layout_matrix=rbind(cbind(1,2),cbind(3,4)),top = textGrob(paste0(sel.season,", SAOD (SVD1 SST-SLP): Obs & EC-Earth3 DCPP-LF/HF"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_EC1_DCPP_OBS_r_LF_HF.png"),plot=fig,width = 13, height = 7)

