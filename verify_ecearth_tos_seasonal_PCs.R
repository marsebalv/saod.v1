# This script compares the PCs  of the initialized decadal predictions of SST (tos) of 
# EC-Earth CMIP6 with observed data
#
# M. Alvarez (2020)
#--------------------------------------------------------------------------------
rm(list=ls())
graphics.off()

setwd("/home/maralv/")

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

#---------------------------------------------------------------------------------------
#  Man Program
#---------------------------------------------------------------------------------------
for(sel.lead in 1:10){
  
#################### Settings ########################
# Remove trend? TRUE or FALSE
remove.trend=TRUE
# Select season
sel.season="JJA"
# # Lead selection
# sel.lead=1
######################################################


#_________________________________________
# Loads data                              \_____________________________________________

# Observations 

if(remove.trend==TRUE){
  load(file=paste0("/home/maralv/data/",sel.season,"_obs_PCs_EOFs_SST_weighted_notrend.RData"))
}else{
  load(file=paste0("/home/maralv/data/",sel.season,"_obs_PCs_EOFs_SST_weighted.RData"))  
} 
sst.pcs$normpc2=NULL
sst.pcs$normpc3=NULL
setnames(sst.pcs,"normpc1","normpc1.obs")

obs.sst.pcs=sst.pcs
obs.sst.eof=sst.eof
rm(sst.pcs,sst.eof)

if(sel.season == "DJF"){
  mmdd="-01-01"
}else if(sel.season == "MAM"){
  mmdd="-03-01"
}else if(sel.season == "JJA"){
  mmdd="-06-01"
}else if(sel.season == "SON"){
  mmdd="-09-01"
}

# EC-Earth3
if(remove.trend==TRUE){
  load(file=paste0("/home/maralv/data/",sel.season,"_lead",sel.lead,"_ECEarth3_PCs_EOFs_SST_weighted_notrend.RData"))
}else{
  load(file=paste0("/home/maralv/data/",sel.season,"_lead",sel.lead,"_ECEarth3_PCs_EOFs_SST_weighted.RData"))
} 

# Change all dates to day=01 to merge in same date
# Add start date as a column
sst.pcs[,date := as.Date(as.character(paste0(year(targetdate),mmdd)))]
sst.pcs$targetdate=NULL

vfy.sst.pcs=merge(obs.sst.pcs,sst.pcs[,c("date","pc1")],by="date",all=TRUE)
setnames(vfy.sst.pcs,"pc1","normpc1.ece.lead")

# Compute temporal correlation between PCs
temp.cor=cor(vfy.sst.pcs$normpc1.obs,vfy.sst.pcs$normpc1.ece.lead,use = "pairwise.complete.obs")
p.value.temp=cor.test(vfy.sst.pcs$normpc1.obs,vfy.sst.pcs$normpc1.ece.lead,method="pearson",alternative="two.sided",conf.level=0.95)$p.value

# Compute spatial correlation between EOFs
aux=merge(obs.sst.eof,sst.eof,by=c("lat","lon"))
spat.cor=cor(aux$hcorrmap.eof1.x,aux$hcorrmap.eof1.y,use = "pairwise.complete.obs")
p.value.spat=cor.test(aux$hcorrmap.eof1.x,aux$hcorrmap.eof1.y,method="pearson",alternative="two.sided",conf.level=0.95)$p.value

# Adjust sign of correlations coefficients and EOFs according to sign of spat cor.
c=1
if(spat.cor<0){
  c=-1
}
temp.cor=c*temp.cor
spat.cor=c*spat.cor
#------

# Plot

fig.pcs <- ggplot() +
  theme_bw()+
  geom_line(data=vfy.sst.pcs,aes(date, normpc1.obs))+
  geom_line(data=vfy.sst.pcs,aes(date, c*normpc1.ece.lead),color="red")+
  scale_x_date(limits=c(vfy.sst.pcs[year(date)==1960]$date,vfy.sst.pcs[year(date)==2018]$date),breaks=vfy.sst.pcs$date[seq(1,length(vfy.sst.pcs$date),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1,size=12))+
  xlab("Date")+ylab("PC1 (normalized units)")+
  theme(axis.title = element_text(size=12),title = element_text(size=11))+
  if(abs(temp.cor)>p.value.temp){
    # correlation is significantly different from 0
    labs(title = paste0("Principal components obs (black) and EC-Earth3 (red), lead = ",sel.lead,". r.pc = ",as.character(round(temp.cor,3)),"*"))
  }else{
    labs(title = paste0("Principal components obs (black) and EC-Earth3 (red), lead = ",sel.lead,". r.pc = ",as.character(round(temp.cor,3))))
  }


# EOFs figures

bmin=-0.6
bmax=0.6
bstep=0.1
bbreaks=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6)
map.world <- map_data ("world2", wrap = c(-180,180))

fig.eof.obs <- ggplot() +
  geom_contour_fill(data=obs.sst.eof,aes(lon, lat, z = hcorrmap.eof1),na.fill=TRUE)+
  scale_fill_distiller(name="EOF1",palette="RdBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0(sel.season," EOF1 SST. Observations."))+
  theme(axis.text=element_text(size=12),title = element_text(size=11))

fig.eof.ece <- ggplot() +
  geom_contour_fill(data=sst.eof,aes(lon, lat, z = c*hcorrmap.eof1),na.fill=TRUE)+
  scale_fill_distiller(name="EOF1",palette="RdBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  theme(axis.text=element_text(size=12),title = element_text(size=11))+
  if(abs(spat.cor)>p.value.spat){
    # correlation is significantly different from 0
    ggtitle(paste0(sel.season," EOF1 SST. EC-Earth3. Lead = ",as.character(sel.lead),". r.eof = ",as.character(round(spat.cor,3)),"*"))
  }else{
    ggtitle(paste0(sel.season, "EOF1 SST. EC-Earth3. Lead = ",as.character(sel.lead),". r.eof = ",as.character(round(spat.cor,3))))
  }
  
fig = grid.arrange(fig.eof.obs,fig.eof.ece,fig.pcs, layout_matrix = cbind(c(1,3), c(2,3)))

if(remove.trend==TRUE){
  ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_verify_EOF_PCs_lead_",sel.lead,"_SST_weighted_notrend.png"),plot=fig,width = 11, height = 8)
}else{
  ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_verify_EOF_PCs_lead_",sel.lead,"SST_weighted.png"),plot=fig,width = 11, height = 8)
} 
rm(list=ls())
}

