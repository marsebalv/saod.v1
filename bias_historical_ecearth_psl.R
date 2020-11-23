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

#_________________________________________
# Opens data files into single data table \_____________________________________________

# Data has been regridded to a regular 2,5x2,5 regular grid using cdo:
# #!/bin/bash
# cd /pathtodata
#
# for pslfile in *
#   do
# echo $pslfile
# cdo remapbil,r144x73 /pathtodata/$pslfile /pathnewdata/${pslfile::-3}_r144x73.nc
# done

ipath="/home/maralv/data/"

# Define subset of lat/lon
subset = list(lat = list(-60:60), lon = list(0:20,250:360),time = c("1949-01-01", "2014-12-31"))
members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)

# Only for the first member
mmb=1
# Form input file name
ifile=paste0(ipath,"psl_Amon_EC-Earth3_historical_r",mmb,"i1p1f1_gr_185001-201412_r144x73.nc")
# Open file as data.frame
data1 = metR::ReadNetCDF(ifile,vars=c("psl"), out='data.frame',subset = subset)
# Modify name of variables
colnames(data1) = c("date","lat","lon",paste0("psl.",mmb))

# For the rest of the members, merge sst with the previous data.frame
for(mmb in members[2:15]){
  ifile=paste0(ipath,"psl_Amon_EC-Earth3_historical_r",mmb,"i1p1f1_gr_185001-201412_r144x73.nc")
  data = metR::ReadNetCDF(ifile,vars=c("psl"), out='data.frame',subset = subset)
  colnames(data) = c("date","lat","lon",paste0("psl.",mmb))
  data1=merge(data1,data,by=c("date","lat","lon"),all=TRUE)
}
rm(data)
as.data.table(data1)

psl.ECE=data1

psl.ECE$date=as.Date(psl.ECE$date)
rm(data1)

# Change Format of dates and day always 01.
psl.ECE$date = as.Date(paste(format(psl.ECE$date,'%Y-%m'),"-01",sep=""))

# Change units to hPa
for(mmb in members){
  eval(parse(text = paste0("psl.ECE$psl.",mmb, "=","psl.ECE$psl.",mmb,"/100" )))
}

#_________________________
# Computes Ensemble Mean  \______________________________________________

# Computes ensemble mean
psl.ECE[, psl.em := mean(as.matrix(.SD), na.rm = TRUE), by = .(lon, lat, date)]

# Add month
psl.ECE[, targetmonth := (month(date))]


#_________________________
# Opens observations      \______________________________________________
# Sea level pressure from NCEP/NCAR

# Jan 1854 - Sep 2016
file.in="/home/maralv/data/slp.ncepncar.mon.nc"

# Longitude regions (divided in the Atlantic by Greenwich)
subset.in = list(Y = list(-60:60), X = list(0:20,250:360), T = c("1949-01-10", "2014-12-31"))


psl = ReadNetCDF(file.in, vars = "pressure", out = "data.frame",
                  subset = subset.in, key = TRUE)


# Rearrange column names, formats
colnames(psl) = c("date","lat","lon","psl.obs")

# Adjust dimensions of date
psl$date = as.Date(psl$date)

# Change Format of dates and day always 01.
psl$date = as.Date(paste(format(psl$date,'%Y-%m'),"-01",sep=""))

# Change to data table
psl=as.data.table(psl)

# Change to hPa
psl$psl.obs=psl$psl.obs/100

#________________________________________
# Merge & save data table                \______________________________________________

psl.all=merge(psl,psl.ECE,by=c("date","lat","lon"),all=TRUE)

# # Apply mask for values over land
# # Define mask
# mask2 <- psl.all[, .(lon = lon, lat = lat, land = MaskLand(lon, lat))]
# # Apply mask
# psl.all$land=mask2$land
# psl.all = psl.all[land==FALSE,]
# 
# psl.all$land = NULL

# save("psl.all",file="/home/maralv/data/full.obs.psl.and.hist.ECEarth3.RData")

#________________________________________
# Load data                              \______________________________________________

# load(file="/home/maralv/data/full.obs.psl.and.hist.ECEarth3.RData")

#################### Settings ########################
# Select season
sel.season="JJA"
# Members being used
members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)
######################################################

#________________________________________
# Perform seasonal averages              \______________________________________________

# Define seasons
psl.all = psl.all[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
psl.all = psl.all[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
psl.all = psl.all[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
psl.all = psl.all[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

psl.all = psl.all[,season.year := year(date)]

# Adjust december to be part of season.year +1 year (year of JF)
psl.all[month(date)==12,]$season.year=psl.all[month(date)==12,]$season.year+1
psl.all = psl.all[season == sel.season,]

# Perform seasonal averages
for (mmb in members) {
  psl.all[ , eval(parse(text = paste0("psl.",mmb, ":=ave(psl.",mmb,")"))),by=.(season.year,lat,lon)]
}
psl.all = psl.all[,psl.em := ave(psl.em),by=.(season.year,lat,lon)]

# Remove repeated rows
psl.all = unique(psl.all, by=c("lat","lon","season.year"))

# Change longitude to -180:180 to plot correctly
psl.all[lon>180]$lon=psl.all[lon>180]$lon-360

#__________________________________________
# Compute bias                             \______________________________________________

# Compute bias for each member
for (mmb in members) {
  psl.all[ , eval(parse(text = paste0("bias.psl.",mmb, ":=","psl.",mmb,"-psl.obs" )))]
}
# Compute bias for the ensemble mean
psl.all[,bias.psl.em:=psl.em-psl.obs]

# Average over all dates
cols = c("bias.psl.1","bias.psl.2","bias.psl.4","bias.psl.5","bias.psl.6","bias.psl.7","bias.psl.8","bias.psl.9","bias.psl.10","bias.psl.11","bias.psl.12","bias.psl.13","bias.psl.14","bias.psl.15","bias.psl.16","bias.psl.em")
psl.ave.bias =  psl.all[, lapply(.SD, mean,na.rm=TRUE), .SDcols=cols, by=.(lat,lon)]

# Now I want the bias in the whole basin for each season, and in the two rectangles NE and SW
psl.season.basin.bias = psl.all[lat>-46 & lat<0 & lon>-50 & lon<20,]
psl.season.N.bias = psl.all[lat>-46 & lat<(-23) & lon>-50 & lon<20,]
psl.season.S.bias = psl.all[lat>(-23) & lat<(0) & lon>-50 & lon<20,]

psl.season.basin.bias =  psl.season.basin.bias[, lapply(.SD, mean,na.rm=TRUE), .SDcols=cols, by=.(season.year)]
psl.season.N.bias =  psl.season.N.bias[, lapply(.SD, mean,na.rm=TRUE), .SDcols=cols, by=.(season.year)]
psl.season.S.bias =  psl.season.S.bias[, lapply(.SD, mean,na.rm=TRUE), .SDcols=cols, by=.(season.year)]

# Compute max and min across all ensemble members biases
psl.season.basin.bias = psl.season.basin.bias[, maxi := max(.SD),.SDcols=cols, by=season.year]
psl.season.basin.bias = psl.season.basin.bias[, mini := min(.SD),.SDcols=cols, by=season.year]

psl.season.N.bias = psl.season.N.bias[, maxi := max(.SD),.SDcols=cols, by=season.year]
psl.season.N.bias = psl.season.N.bias[, mini := min(.SD),.SDcols=cols, by=season.year]

psl.season.S.bias = psl.season.S.bias[, maxi := max(.SD),.SDcols=cols, by=season.year]
psl.season.S.bias = psl.season.S.bias[, mini := min(.SD),.SDcols=cols, by=season.year]

#__________________________________________
# Plots and maps                           \______________________________________________
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-3.5
bmax=3.5
bstep=0.5
bbreaks=seq(bmin,bmax,bstep)

g1<- ggplot(data=psl.ave.bias) +
  geom_contour_fill(aes(lon, lat, z = bias.psl.em),breaks=bbreaks,na.fill=TRUE)+
  scale_fill_distiller(name="bias",palette="RdYlBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-100,20,20))+
  scale_y_latitude(breaks=seq(-60,60,20))+
  geom_rect(aes(xmin=-50, xmax=20, ymin=-46, ymax=0),fill=NA, color="black", alpha=0.05,linetype="dashed")+
  geom_rect(aes(xmin=-49.5, xmax=19.5, ymin=-45.5, ymax=-23.25),fill=NA, color="#4dac26", alpha=0.05,linetype="solid",size=0.25)+
  geom_rect(aes(xmin=-49.5, xmax=19.5, ymin=-22.75, ymax=-0.5),fill=NA, color="#5e3c99", alpha=0.05,linetype="solid",size=0.25)+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="black", color="black", inherit.aes = F)+
  ggtitle(paste0("EC-Earth3 historical ensemble mean bias. Season ",sel.season))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 <- ggplot(psl.season.basin.bias) +
  theme_bw()+
  geom_line(aes(season.year, bias.psl.em),color="#e66101")+
  geom_ribbon(aes(x=season.year, ymin=mini , ymax=maxi ),fill="#fdb863",alpha=0.4)+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,0.50),limits=c(-3.5,3.5),expand = c(0., 0.))+
  scale_x_continuous(breaks = seq(1950,2015,5),limits=c(1949,2014),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Year of season")+ylab("SLP bias (hPa)")+
  ggtitle(paste0("Average over the South Atlantic Basin (dashed rectangle)"))+
  theme(axis.title = element_text(size=12),title = element_text(size=10))



g3 <- ggplot() +
  theme_bw()+
  geom_line(data=psl.season.N.bias,aes(season.year, bias.psl.em),color="#5e3c99")+
  geom_ribbon(data=psl.season.N.bias,aes(x=season.year, ymin=mini , ymax=maxi ),fill="#b2abd2",alpha=0.4)+

  geom_line(data=psl.season.S.bias,aes(season.year, bias.psl.em),color="#4dac26")+
  geom_ribbon(data=psl.season.S.bias,aes(x=season.year, ymin=mini , ymax=maxi ),fill="#b8e186",alpha=0.2)+

  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-6,6,1),limits=c(-6,6),expand = c(0., 0.))+
  scale_x_continuous(breaks = seq(1950,2015,5),limits=c(1949,2014),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Year of season")+ylab("SLP bias (hPa)")+
  ggtitle(paste0("Average over the N (purple) and S (green) regions of South Atlantic Basin"))+
  theme(axis.title = element_text(size=12),title = element_text(size=10))

fig <- grid.arrange(g1,g2,g3,layout_matrix=rbind(c(1,2),c(1,3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_PSL_biases.png"),plot=fig,width = 14, height = 8)

