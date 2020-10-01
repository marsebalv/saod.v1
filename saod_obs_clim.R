# This script opens SST and SLP data in the South Atlantic region and computes
# and plots the climatology for each month
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
#Datos de SLP (NCEP-NCAR Reana1)
# lat e/2.5°
# lon e/2.5°

# File to load
file.in="/home/maralv/data/slp.ncepncar.mon.nc"

# Longitude regions (divided in the Atlantic by Greenwich)
lon1=seq(0,20,2.5)
lon2=seq(290,357.5,2.5)

subset.in1 = list(Y = seq(-40,0,2.5), X=lon1)
subset.in2 = list(Y = seq(-40,0,2.5), X=lon2)
# Even though data is listed from 0 to -40 it makes no difference (at least reading as a data frame)

# Load data in data.frames
slp1 = ReadNetCDF(file.in, vars = "pressure", out = "data.frame",
                 subset = subset.in1, key = TRUE)

slp2 = ReadNetCDF(file.in, vars = "pressure", out = "data.frame",
                  subset = subset.in2, key = TRUE)

# Combine regions in one dataframe
slp=rbind(slp1,slp2)
rm(slp1,slp2)

colnames(slp) = c("date","lat","lon","slp")


# Change slp units to hPa
slp$slp = slp$slp/100

#----------

# NOAA Extended Reconstructed Sea Surface Temperature (SST) V5

# Jan 1854 - Sep 2016
file.in="/home/maralv/data/sst.noaaersstv5.mon.nc"

# Longitude regions (divided in the Atlantic by Greenwich)
lon1=seq(0,20,2)
lon2=seq(290,358,2)

subset.in1 = list(Y = seq(-40,0,2), X=lon1)
subset.in2 = list(Y = seq(-40,0,2), X=lon2)

sst1 = ReadNetCDF(file.in, vars = "sst", out = "data.frame",
                  subset = subset.in1, key = TRUE)
sst2 = ReadNetCDF(file.in, vars = "sst", out = "data.frame",
                  subset = subset.in2, key = TRUE)

# Combine regions in one dataframe
sst=rbind(sst1,sst2)
rm(sst1,sst2)

# Rearrange column names, formats
sst$zlev = NULL
colnames(sst) = c("date","lat","lon","sst")

# Adjust dimensions of date
sst$date = as.Date(sst$date)
slp$date = as.Date(slp$date)

# Define years
sst = sst[year(date)>=1949 & year(date)<2020,]
slp = slp[year(date)>=1949 & year(date)<2020,]

# Change Format of dates and day always 01.
sst$date = as.Date(paste(format(sst$date,'%Y-%m'),"-01",sep=""))
# Change Format of dates
slp$date = as.Date(paste(format(slp$date,'%Y-%m'),"-01",sep=""))

# Change to data table
sst=as.data.table(sst)
slp=as.data.table(slp)

# Apply mask for values over land
# Define mask
mask <- slp[, .(lon = lon, lat = lat, land = MaskLand(lon, lat))]
mask2 <- sst[, .(lon = lon, lat = lat, land = MaskLand(lon, lat))]
# Apply mask
slp$land=mask$land
slp = slp[land==FALSE,]
sst$land=mask2$land
sst = sst[land==FALSE,]

sst$land = NULL
slp$land = NULL

#---------------------------------------------------------------------------------------
#  Compute climatology and anomalies 
#---------------------------------------------------------------------------------------

slp$month = month(as.Date(slp$date))
sst$month = month(as.Date(sst$date))
# Average over all years
clim.slp = slp[, .(clim = ave(slp)), by=.(lat,lon,month)]
clim.sst = sst[, .(clim = ave(sst)), by=.(lat,lon,month)]
# Remove repeated rows
clim.slp = unique(clim.slp, by=c("lat","lon","month"))
clim.sst = unique(clim.sst, by=c("lat","lon","month"))

# Prepare variable
clim.slp$runmean3=NA_real_
clim.sst$runmean3=NA_real_

# Ojo, no tengo que pasarle como argumento a la funci+on "clim$clim1", y tengo que repetir,
# más allá del "by", qué otras variables quiero seguir teniendo en la DT
clim.slp = clim.slp[,.(month,smthclim = fun3mrm(clim)),by=.(lat,lon)]

clim.sst = clim.sst[,.(month,smthclim = fun3mrm(clim)),by=.(lat,lon)]

################################################
# Plots climatology according to month

# Change longitude to -180:180 to plot correctly
clim.sst[lon>180]$lon=clim.sst[lon>180]$lon-360
clim.slp[lon>180]$lon=clim.slp[lon>180]$lon-360

# Plot a sample
map.world <- map_data ("world2", wrap = c(-180,180))

# Make month a factor to reorder facets
clim.sst$month_f = factor(clim.sst$month, levels=c('12','1','2','3','4','5','6','7','8','9','10','11'))
clim.slp$month_f = factor(clim.slp$month, levels=c('12','1','2','3','4','5','6','7','8','9','10','11'))

month.labs <- c("December","January","February","March","April","May","June","July","August","September","October","November")
names(month.labs) <- c('12','1','2','3','4','5','6','7','8','9','10','11')

# Sin ponerle color a los contornos
ggplot() +
  geom_contour_fill(data=clim.sst,aes(lon, lat, z = smthclim),breaks=seq(10,30,2.5),na.fill=TRUE)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       breaks=seq(10,30,2.5),
                       limits=c(10,30),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  new_scale_color() +
  geom_contour(data=clim.slp,aes(lon, lat, z = smthclim),breaks=seq(1010,1024,2),color="black",size=0.25)+
  scale_color_distiller(name="SLP",palette="Greys",direction=1,
                        breaks=seq(1010,1024,2),
                        limits=c(1010,1024))+
  geom_text_contour(data=clim.slp,aes(lon, lat, z = smthclim),breaks=seq(1010,1024,2),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  theme(axis.text=element_text(size=12))+
  facet_wrap(~ month_f, ncol=3,labeller = labeller(month_f = month.labs)) +
  labs(title="SLP (NCEP/NCAR) and SST (ERSSTv5) 1949-2019 monthly climatology")+
  theme(strip.background = element_rect(color="black", fill="white", size=0.6, linetype="solid"))+
  theme(strip.text = element_text(size = 12, colour = "black"))

ggsave("/home/maralv/Dropbox/DMI/Figures/clima_monthly_SLP_SST_4919.png",plot=last_plot(),width = 10, height = 8)

#########################

rm(mask,mask2,subset.in1,subset.in2,file.in,lon1,lon2)
