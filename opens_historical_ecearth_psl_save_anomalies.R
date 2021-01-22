# This script opens the initialized decadal predictions of SLP (psl) of 
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
subset = list(lat = list(-45:0), lon = list(0:20,310:357.5),time = c("1949-01-01", "2014-12-31"))
members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)
  
# Only for the first member
mmb=1
# Form input file name
ifile=paste0(ipath,"psl_Amon_EC-Earth3_historical_r",mmb,"i1p1f1_gr_185001-201412_r144x73.nc")
# Open file as data.frame
data1 = metR::ReadNetCDF(ifile,vars=c("psl"), out='data.frame',subset = subset)
# Modify name of variables
colnames(data1) = c("targetdate","lat","lon",paste0("psl.",mmb))

# For the rest of the members, merge sst with the previous data.frame
for(mmb in members[2:15]){
  ifile=paste0(ipath,"psl_Amon_EC-Earth3_historical_r",mmb,"i1p1f1_gr_185001-201412_r144x73.nc")
  data = metR::ReadNetCDF(ifile,vars=c("psl"), out='data.frame',subset = subset)
  colnames(data) = c("targetdate","lat","lon",paste0("psl.",mmb))
  data1=merge(data1,data,by=c("targetdate","lat","lon"),all=TRUE)
}
rm(data)
as.data.table(data1)

psl.ECE=data1

psl.ECE$targetdate=as.Date(psl.ECE$targetdate)
rm(data1)

#________________________________________
# Computes Ensemble Mean and Climatology \______________________________________________

# Computes ensemble mean
psl.ECE[, psl.em := mean(as.matrix(.SD), na.rm = TRUE), by = .(lon, lat, targetdate)]

# Add month
psl.ECE[, targetmonth := (month(targetdate))]

# Compute climatology
clim.psl.ECE = psl.ECE[, .(clim = ave(psl.em)), by=.(lat,lon,targetmonth)]
# Remove repeated rows
clim.psl.ECE = unique(clim.psl.ECE, by=c("lat","lon","targetmonth"))
# Prepare variable
clim.psl.ECE$smthclim=NA_real_
# Compute 3 point running mean
clim.psl.ECE = clim.psl.ECE[,.(targetmonth,clim,smthclim = fun3mrm(clim)),by=.(lat,lon)]
clim.psl.ECE$clim=NULL



#________________________________________
# Compute anomalies                      \______________________________________________

# Merge climatology
psl.ECE=merge(psl.ECE,clim.psl.ECE,by=c("lat","lon","targetmonth"),all.x=TRUE)
# Move columns to the left
setcolorder(psl.ECE,neworder=c("lat","lon","targetdate","targetmonth"))

# Compute anomalies for each member
for (mmb in members) {
  psl.ECE[ , eval(parse(text = paste0("apsl.",mmb, ":=","psl.",mmb,"-smthclim" )))]
}
# Compute anomaly for the ensemble mean
psl.ECE[,apsl.em:=psl.em-smthclim]

# Remove raw values from data table
for (mmb in members) {
  eval(parse(text = paste0("psl.ECE$psl.",mmb, "=","NULL" )))
}
psl.ECE$psl.em=NULL

psl.ECE$smthclim=NULL

#________________________________________
# Save data table                        \______________________________________________

psl.ECE=psl.ECE[order(targetdate),]

save("psl.ECE",file="/home/maralv/data/apsl.ECEarth.hist.19492014.RData")

save("clim.psl.ECE",file="/home/maralv/data/clim.psl.ECEarth.hist.19492014.RData")
