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


# Create data table to save correlation coefficients
# SST
corr.dcpp.sst = data.table(member=1:16)
corr.dcpp.sst[member==16]$member=90
corr.dcpp.sst = corr.dcpp.sst[, .(lead = 1:10), by = c(colnames(corr.dcpp.sst))]
corr.dcpp.sst$r = NA_real_
corr.hist.sst = data.table(member=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,90))
corr.hist.sst$r = NA_real_
# SLP
corr.dcpp.slp = data.table(member=1:16)
corr.dcpp.slp[member==16]$member=90
corr.dcpp.slp = corr.dcpp.slp[, .(lead = 1:10), by = c(colnames(corr.dcpp.slp))]
corr.dcpp.slp$r = NA_real_
corr.hist.slp = data.table(member=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,90))
corr.hist.slp$r = NA_real_

##############
# Start for on leads
for(l in 1:10){

sel.lead=l  

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

# Plot
breaks.dates=exp.coef.norm.OBS$date
limits.EC=4

g4 <- ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.separated.dcpp[member==1,],aes(targetdate, sst.med),col="#de77ae",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated.dcpp[member==1,], aes(targetdate, ymin=sst.min , ymax=sst.max ),fill="#de77ae",alpha=0.2)+
  geom_ribbon(data=exp.coef.separated.dcpp[member==1,], aes(targetdate, ymin=sst.p20 , ymax=sst.p80 ),fill="#de77ae",alpha=0.4)+
  
  geom_line(data=exp.coef.norm.OBS,aes(date, ec.sst.1),col="black",alpha=0.8,size=0.7)+
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(ceiling(-limits.EC),floor(limits.EC),1),limits=c(-limits.EC,limits.EC),expand = c(0, 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("targetdate")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1 DCPP lead ",sel.lead, ", SST (pink). Dark shading: 9 less extreme mmbs"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

g5 <- ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.separated.dcpp[member==1,],aes(targetdate, slp.med),col="#7fbc41",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated.dcpp[member==1,], aes(targetdate, ymin=slp.min , ymax=slp.max ),fill="#7fbc41",alpha=0.2)+
  geom_ribbon(data=exp.coef.separated.dcpp[member==1,], aes(targetdate, ymin=slp.p20 , ymax=slp.p80 ),fill="#7fbc41",alpha=0.4)+
  
  geom_line(data=exp.coef.norm.OBS,aes(date, ec.slp.1),col="black",alpha=0.8,size=0.7)+
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(ceiling(-limits.EC),floor(limits.EC),1),limits=c(-limits.EC,limits.EC),expand = c(0, 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("targetdate")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1 DCPP lead ",sel.lead, ", SLP (green). Dark shading: 9 less extreme mmbs"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

g6 <- ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.separated.hist[member==1,],aes(date, sst.med),col="#de77ae",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated.hist[member==1,], aes(date, ymin=sst.min , ymax=sst.max ),fill="#de77ae",alpha=0.2)+
  geom_ribbon(data=exp.coef.separated.hist[member==1,], aes(date, ymin=sst.p20 , ymax=sst.p80 ),fill="#de77ae",alpha=0.4)+
  
  geom_line(data=exp.coef.norm.OBS,aes(date, ec.sst.1),col="black",alpha=0.8,size=0.7)+
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(ceiling(-limits.EC),floor(limits.EC),1),limits=c(-limits.EC,limits.EC),expand = c(0, 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1 Historical SST (pink). Dark shading: 9 less extreme mmbs"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

g7 <- ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.separated.hist[member==1,],aes(date, slp.med),col="#7fbc41",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated.hist[member==1,], aes(date, ymin=slp.min , ymax=slp.max ),fill="#7fbc41",alpha=0.2)+
  geom_ribbon(data=exp.coef.separated.hist[member==1,], aes(date, ymin=slp.p20 , ymax=slp.p80 ),fill="#7fbc41",alpha=0.4)+
  
  geom_line(data=exp.coef.norm.OBS,aes(date, ec.slp.1),col="black",alpha=0.8,size=0.7)+
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(ceiling(-limits.EC),floor(limits.EC),1),limits=c(-limits.EC,limits.EC),expand = c(0, 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1 Historical SLP (green). Dark shading: 9 less extreme mmbs"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

#### Plot of spatial fields
map.world <- map_data ("world2", wrap = c(-180,180))

# Only one plot, both patterns
bmin=-0.9
bmax=0.9
bbreaks=c(-99,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,99)

g1 <- ggplot() +
  geom_contour_fill(data=sst.OBS[date==(sst.OBS$date[1])],aes(lon, lat, z = hcm.1, fill=stat(level)),breaks=bbreaks,na.fill=TRUE)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       limits=c(bmin,bmax),
                       super = ScaleDiscretised,
                       guide = guide_colorsteps(),
                       oob  = scales::squish)+
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 10))+
  new_scale_color() +
  geom_contour(data=psl.OBS[date==(sst.OBS$date[1])],aes(lon, lat, z = hcm.1),bbreaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=psl.OBS[date==(psl.OBS$date[1])],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("Observations"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 <- ggplot() +
  geom_contour_fill(data=sst.ECE.allin1.hist[date==sst.ECE.allin1.hist$date[1] & member==1],aes(lon, lat, z = hcm.1, fill=stat(level)),breaks=bbreaks,na.fill=TRUE)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       limits=c(bmin,bmax),
                       super = ScaleDiscretised,
                       guide = guide_colorsteps(),
                       oob  = scales::squish)+
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 10))+
  new_scale_color() +
  geom_contour(data=psl.ECE.allin1.hist[date==psl.ECE.allin1.hist$date[1] & member==1],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=psl.ECE.allin1.hist[date==psl.ECE.allin1.hist$date[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("Historical"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g3 <- ggplot() +
  geom_contour_fill(data=sst.ECE.allin1.dcpp[targetdate==sst.ECE.allin1.dcpp$targetdate[1] & member==1],aes(lon, lat, z = hcm.1, fill=stat(level)),breaks=bbreaks,na.fill=TRUE)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       limits=c(bmin,bmax),
                       super = ScaleDiscretised,
                       guide = guide_colorsteps(),
                       oob  = scales::squish)+
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 10))+
  new_scale_color() +
  geom_contour(data=psl.ECE.allin1.dcpp[targetdate==psl.ECE.allin1.dcpp$targetdate[1] & member==1],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=psl.ECE.allin1.dcpp[targetdate==psl.ECE.allin1.dcpp$targetdate[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("Decadal Prediction. Lead: ",sel.lead))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

# Add in 1 figure
lmatx=rbind(c(1,1,2,2,3,3),cbind(4,4,4,5,5,5),cbind(6,6,6,7,7,7))

fig <- grid.arrange(g1,g2,g3,g6,g4,g7,g5, layout_matrix=lmatx,top = textGrob(paste0(sel.season,", SAOD (SVD1 SST-SLP): Obs & EC-Earth3"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_SAOD_obs-hist-dcpp_allmembers_lead_year_",sel.lead,".png"),plot=fig,width = 12, height = 8)

##########################################
# Compute correlation between obs & dcpp

setnames(exp.coef.separated.dcpp,"targetdate","date")

a=merge(exp.coef.norm.OBS[,.(date,ec.sst.1,ec.slp.1)],exp.coef.separated.dcpp[,.(member,date,ec1.sst.norm,ec1.slp.norm,sst.med,slp.med)],by=c("date"))
setnames(a,c("ec1.sst.norm","ec1.slp.norm","sst.med","slp.med"),c("ec1.sst.norm.dcpp","ec1.slp.norm.dcpp","sst.med.dcpp","slp.med.dcpp"))
aux=a[,.(date,member,ec.sst.1,ec1.sst.norm.dcpp)]
aux= aux[, r := cor(ec.sst.1,ec1.sst.norm.dcpp,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
# Save correlation coefficients
corr.dcpp.sst[lead==sel.lead & member<=15,]$r=aux$r
# With the mean EC1
aux=unique(a[,.(date,ec.sst.1,sst.med.dcpp)])
corr.dcpp.sst[lead==sel.lead & member==90,]$r=cor(aux$ec.sst.1,aux$sst.med.dcpp,use="pairwise.complete.obs")

rm(aux)
aux=a[,.(date,member,ec.slp.1,ec1.slp.norm.dcpp)]
aux= aux[, r := cor(ec.slp.1,ec1.slp.norm.dcpp,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
# Save correlation coefficients
corr.dcpp.slp[lead==sel.lead & member<=15,]$r=aux$r
# With the mean EC1
aux=unique(a[,.(date,ec.slp.1,slp.med.dcpp)])
corr.dcpp.slp[lead==sel.lead & member==90,]$r=cor(aux$ec.slp.1,aux$slp.med.dcpp,use="pairwise.complete.obs")


# Remove dcpp data
rm(sst.ECE.allin1.dcpp,psl.ECE.allin1.dcpp,a)
} # End for on leads

##########################################
# Compute correlation between obs & hist

# Create a dt merging
a=merge(exp.coef.norm.OBS[,.(date,ec.sst.1,ec.slp.1)],exp.coef.separated.hist[,.(member,date,ec1.sst.norm,ec1.slp.norm,sst.med,slp.med)],by="date")
setnames(a,c("ec1.sst.norm","ec1.slp.norm","sst.med","slp.med"),c("ec1.sst.norm.hist","ec1.slp.norm.hist","sst.med.hist","slp.med.hist"))
# Reduce to selection for sst
aux=a[,.(date,member,ec.sst.1,ec1.sst.norm.hist)]
aux= aux[, r := cor(ec.sst.1,ec1.sst.norm.hist,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
# Save correlation coefficients
corr.hist.sst[1:15]$r=aux$r

aux=unique(a[,.(date,ec.sst.1,sst.med.hist)])
corr.hist.sst[16]$r=cor(aux$ec.sst.1,aux$sst.med.hist,use="pairwise.complete.obs")
rm(aux)
# Reduce to selection for slp
aux=a[,.(date,member,ec.slp.1,ec1.slp.norm.hist)]
aux= aux[, r := cor(ec.slp.1,ec1.slp.norm.hist,use="pairwise.complete.obs"),by="member"]
aux = unique(aux[,.(member,r)])
# Save correlation coefficients
corr.hist.slp[1:15]$r=aux$r

aux=unique(a[,.(date,ec.slp.1,slp.med.hist)])
corr.hist.slp[16]$r=cor(aux$ec.slp.1,aux$slp.med.hist,use="pairwise.complete.obs")



##################
# Plotting
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

h1 <- ggplot()+
  theme_bw()+
  geom_point(data=corr.dcpp.sst[member<=16,],aes(x=lead,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_line(data=corr.dcpp.sst[member==90,],aes(x=lead,y=r))+
  scale_x_continuous(breaks=seq(1,10,1),expand = c(0.02, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-0.4,0.4,0.1),limits=c(-0.45,0.45),expand = c(0, 0.))+
  xlab("Lead")+ylab("Correlation")+
  ggtitle(paste0("Correlation between OBS & DCPP EC1-SST"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))

h2 <-ggplot()+
  theme_bw()+
  geom_point(data=corr.hist.sst[member<=16,],aes(x=1,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_point(data=corr.hist.sst[member==90,],aes(x=1,y=r),shape=3,size=10)+
  scale_x_continuous(breaks=1,expand = c(0.02, 0),labels=c("Historical"))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-0.4,0.4,0.1),limits=c(-0.45,0.45),expand = c(0, 0.))+
  ylab("Correlation")+xlab("All")+
  ggtitle(paste0("& HIST"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))


h3 <- ggplot()+
  theme_bw()+
  geom_point(data=corr.dcpp.slp[member<=16,],aes(x=lead,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_line(data=corr.dcpp.slp[member==90,],aes(x=lead,y=r))+
  scale_x_continuous(breaks=seq(1,10,1),expand = c(0.02, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-0.4,0.4,0.1),limits=c(-0.45,0.45),expand = c(0, 0.))+
  xlab("Lead")+ylab("Correlation")+
  ggtitle(paste0("Correlation between OBS & DCPP EC1-SLP"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))

h4 <-ggplot()+
  theme_bw()+
  geom_point(data=corr.hist.slp[member<=16,],aes(x=1,y=r,group=member,color=as.factor(member)))+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")+
  geom_point(data=corr.hist.slp[member==90,],aes(x=1,y=r),shape=3,size=10)+
  scale_x_continuous(breaks=1,expand = c(0.02, 0),labels=c("Historical"))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-0.4,0.4,0.1),limits=c(-0.45,0.45),expand = c(0, 0.))+
  ylab("Correlation")+xlab("All")+
  ggtitle(paste0("& HIST"))+
  theme(axis.title = element_text(size=13),title = element_text(size=12))


fig <- grid.arrange(h1,h2,h3,h4, layout_matrix=rbind(cbind(1,1,1,1,2),cbind(3,3,3,3,4)),top = textGrob(paste0(sel.season,", SAOD (SVD1 SST-SLP): Obs & EC-Earth3"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_EC1_DCPP_HIST_OBS_r.png"),plot=fig,width = 9, height = 6)

