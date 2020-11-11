# This script opens the initialized decadal predictions of SST (tos) of 
# EC-Earth CMIP6
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
#  functions 
#---------------------------------------------------------------------------------------

# Apply a 3-month moving running mean to monthly data to smooth climatology
fun3mrm <- function(datos) {
  # Repeat month for recursiveness
  a = c(datos[length(datos)],datos,datos[1])
  rdo = rollmean(a,3,align="center")
  return(rdo)
}

#---------------------------------------------------------------------------------------
#  Man Program
#---------------------------------------------------------------------------------------

# #_________________________________________
# # Opens data files into single data table \_____________________________________________
# 
# # Data has been regridded to a regular 2,5x2,5 regular grid using cdo:
# # #!/bin/bash
# # cd /pathtodata
# # 
# # for pslfile in *
# #   do
# # echo $pslfile
# # cdo remapbil,r144x73 /pathtodata/$pslfile /pathnewdata/${pslfile::-3}_r144x73.nc
# # done
# 
# ipath="/home/maralv/data/"
# 
# # Define subset of lat/lon
# subset = list(lat = list(-60:60), lon = list(0:20,250:360),time = c("1949-01-10", "2014-12-31"))
# members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)
#   
# # Only for the first member
# mmb=1
# # Form input file name
# ifile=paste0(ipath,"tos_Omon_EC-Earth3_historical_r",mmb,"i1p1f1_gn_185001-201412_r180x91.nc")
# # Open file as data.frame
# data1 = metR::ReadNetCDF(ifile,vars=c("tos"), out='data.frame',subset = subset)
# # Modify name of variables
# colnames(data1) = c("date","lat","lon",paste0("sst.",mmb))
# 
# # For the rest of the members, merge sst with the previous data.frame
# for(mmb in members[2:15]){
#   ifile=paste0(ipath,"tos_Omon_EC-Earth3_historical_r",mmb,"i1p1f1_gn_185001-201412_r180x91.nc")
#   data = metR::ReadNetCDF(ifile,vars=c("tos"), out='data.frame',subset = subset)
#   colnames(data) = c("date","lat","lon",paste0("sst.",mmb))
#   data1=merge(data1,data,by=c("date","lat","lon"),all=TRUE)
# }
# rm(data)
# as.data.table(data1)
# 
# sst.ECE=data1
# 
# sst.ECE$date=as.Date(sst.ECE$date)
# rm(data1)
# 
# # Change Format of dates and day always 01.
# sst.ECE$date = as.Date(paste(format(sst.ECE$date,'%Y-%m'),"-01",sep=""))
# 
# #_________________________
# # Computes Ensemble Mean  \______________________________________________
# 
# # Computes ensemble mean
# sst.ECE[, sst.em := mean(as.matrix(.SD), na.rm = TRUE), by = .(lon, lat, date)]
# 
# # Add month
# sst.ECE[, targetmonth := (month(date))]
# 
# 
# #_________________________
# # Opens observations      \______________________________________________
# # NOAA Extended Reconstructed Sea Surface Temperature (SST) V5
# 
# # Jan 1854 - Sep 2016
# file.in="/home/maralv/data/sst.noaaersstv5.mon.nc"
# 
# # Longitude regions (divided in the Atlantic by Greenwich)
# subset.in = list(Y = list(-60:60), X = list(0:20,250:360), T = c("1949-01-10", "2014-12-31"))
# 
# 
# sst = ReadNetCDF(file.in, vars = "sst", out = "data.frame",
#                   subset = subset.in, key = TRUE)
# 
# 
# # Rearrange column names, formats
# sst$zlev = NULL
# colnames(sst) = c("date","lat","lon","sst.obs")
# 
# # Adjust dimensions of date
# sst$date = as.Date(sst$date)
# 
# # Change Format of dates and day always 01.
# sst$date = as.Date(paste(format(sst$date,'%Y-%m'),"-01",sep=""))
# 
# # Change to data table
# sst=as.data.table(sst)
# 
# #________________________________________
# # Merge & save data table                \______________________________________________
# 
# sst.all=merge(sst,sst.ECE,by=c("date","lat","lon"),all=TRUE)
# 
# # Apply mask for values over land
# # Define mask
# mask2 <- sst.all[, .(lon = lon, lat = lat, land = MaskLand(lon, lat))]
# # Apply mask
# sst.all$land=mask2$land
# sst.all = sst.all[land==FALSE,]
# 
# sst.all$land = NULL
# 
# save("sst.all",file="/home/maralv/data/full.obs.sst.and.hist.ECEarth3.RData")

#________________________________________
# Load data                              \______________________________________________

load(file="/home/maralv/data/full.obs.sst.and.hist.ECEarth3.RData")

#################### Settings ########################
# Select season
sel.season="DJF"
# Members being used
members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)
######################################################

#________________________________________
# Perform seasonal averages              \______________________________________________

# Define seasons
sst.all = sst.all[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
sst.all = sst.all[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
sst.all = sst.all[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
sst.all = sst.all[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

sst.all = sst.all[,season.year := year(date)]

# Adjust december to be part of season.year +1 year (year of JF)
sst.all[month(date)==12,]$season.year=sst.all[month(date)==12,]$season.year+1
sst.all = sst.all[season == sel.season,]

# Perform seasonal averages
for (mmb in members) {
  sst.all[ , eval(parse(text = paste0("sst.",mmb, ":=ave(sst.",mmb,")"))),by=.(season.year,lat,lon)]
}
sst.all = sst.all[,sst.em := ave(sst.em),by=.(season.year,lat,lon)]

# Remove repeated rows
sst.all = unique(sst.all, by=c("lat","lon","season.year"))

# Change longitude to -180:180 to plot correctly
sst.all[lon>180]$lon=sst.all[lon>180]$lon-360

#__________________________________________
# Compute bias                             \______________________________________________

# Compute bias for each member
for (mmb in members) {
  sst.all[ , eval(parse(text = paste0("bias.sst.",mmb, ":=","sst.",mmb,"-sst.obs" )))]
}
# Compute bias for the ensemble mean
sst.all[,bias.sst.em:=sst.em-sst.obs]

# Average over all dates
cols = c("bias.sst.1","bias.sst.2","bias.sst.4","bias.sst.5","bias.sst.6","bias.sst.7","bias.sst.8","bias.sst.9","bias.sst.10","bias.sst.11","bias.sst.12","bias.sst.13","bias.sst.14","bias.sst.15","bias.sst.16","bias.sst.em")
sst.ave.bias =  sst.all[, lapply(.SD, mean), .SDcols=cols, by=.(lat,lon)]

# Now I want the bias in the whole basin for each season, and in the two rectangles NE and SW
sst.season.basin.bias = sst.all[lat>-46 & lat<0 & lon>-50 & lon<20,]
sst.season.NE.bias = sst.all[lat>-14 & lat<0 & lon>-10 & lon<10,]
sst.season.SW.bias = sst.all[lat>-44 & lat<(-30) & lon>-35 & lon<(-15),]

sst.season.basin.bias =  sst.season.basin.bias[, lapply(.SD, mean,na.rm=TRUE), .SDcols=cols, by=.(season.year)]
sst.season.NE.bias =  sst.season.NE.bias[, lapply(.SD, mean,na.rm=TRUE), .SDcols=cols, by=.(season.year)]
sst.season.SW.bias =  sst.season.SW.bias[, lapply(.SD, mean,na.rm=TRUE), .SDcols=cols, by=.(season.year)]

# Compute max and min across all ensemble members biases
sst.season.basin.bias = sst.season.basin.bias[, maxi := max(.SD),.SDcols=cols, by=season.year]
sst.season.basin.bias = sst.season.basin.bias[, mini := min(.SD),.SDcols=cols, by=season.year]

sst.season.NE.bias = sst.season.NE.bias[, maxi := max(.SD),.SDcols=cols, by=season.year]
sst.season.NE.bias = sst.season.NE.bias[, mini := min(.SD),.SDcols=cols, by=season.year]

sst.season.SW.bias = sst.season.SW.bias[, maxi := max(.SD),.SDcols=cols, by=season.year]
sst.season.SW.bias = sst.season.SW.bias[, mini := min(.SD),.SDcols=cols, by=season.year]

#__________________________________________
# Plots and maps                           \______________________________________________
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-5
bmax=5
bstep=1
bbreaks=seq(bmin,bmax,bstep)
  
g1<- ggplot(data=sst.ave.bias) +
  geom_contour_fill(aes(lon, lat, z = bias.sst.em),breaks=bbreaks,na.fill=TRUE)+
  scale_fill_distiller(name="bias",palette="RdBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-100,20,20))+
  scale_y_latitude(breaks=seq(-60,60,20))+
  geom_rect(aes(xmin=-50, xmax=20, ymin=-46, ymax=0),fill=NA, color="black", alpha=0.05,linetype="dashed")+
  geom_rect(aes(xmin=-35, xmax=-15, ymin=-44, ymax=-30),fill=NA, color="#4dac26", alpha=0.05,linetype="solid",size=0.25)+
  geom_rect(aes(xmin=-10, xmax=10, ymin=-15, ymax=0),fill=NA, color="#5e3c99", alpha=0.05,linetype="solid",size=0.25)+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="black", color="black", inherit.aes = F)+
  ggtitle(paste0("EC-Earth3 historical ensemble mean bias. Season ",sel.season))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 <- ggplot(sst.season.basin.bias) +
  theme_bw()+
  geom_line(aes(season.year, bias.sst.em),color="#e66101")+
  geom_ribbon(aes(x=season.year, ymin=mini , ymax=maxi ),fill="#fdb863",alpha=0.4)+ 
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-1.5,5.5,1),limits=c(-1.5,5.5),expand = c(0., 0.))+
  scale_x_continuous(breaks = seq(1950,2015,5),limits=c(1949,2014),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Year of season")+ylab("SST bias (degrees)")+
  ggtitle(paste0("Average over the South Atlantic Basin (dashed rectangle)"))+
  theme(axis.title = element_text(size=12),title = element_text(size=10))

  

g3 <- ggplot() +
  theme_bw()+
  geom_line(data=sst.season.NE.bias,aes(season.year, bias.sst.em),color="#5e3c99")+
  geom_ribbon(data=sst.season.NE.bias,aes(x=season.year, ymin=mini , ymax=maxi ),fill="#b2abd2",alpha=0.4)+ 
  
  geom_line(data=sst.season.SW.bias,aes(season.year, bias.sst.em),color="#4dac26")+
  geom_ribbon(data=sst.season.SW.bias,aes(x=season.year, ymin=mini , ymax=maxi ),fill="#b8e186",alpha=0.2)+   
  
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-1.5,5.5,1),limits=c(-1.5,5.5),expand = c(0., 0.))+
  scale_x_continuous(breaks = seq(1950,2015,5),limits=c(1949,2014),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Year of season")+ylab("SST bias (degrees)")+
  ggtitle(paste0("Average over the NE (purple) and SW (green) regions of South Atlantic Basin"))+
  theme(axis.title = element_text(size=12),title = element_text(size=10))

fig <- grid.arrange(g1,g2,g3,layout_matrix=rbind(c(1,2),c(1,3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_SST_biases.png"),plot=fig,width = 14, height = 8)

  