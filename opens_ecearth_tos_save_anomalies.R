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

# Help functions of metR
.tidy2matrix <- function(data, formula, value.var, fill = NULL, ...) {
  row.vars <- all.vars(formula[[2]])
  col.vars <- all.vars(formula[[3]])
  data <- data.table::as.data.table(data)
  data[, row__ := .GRP, by = c(row.vars)]
  data[, col__ := .GRP, by = c(col.vars)]
  if (is.null(fill)){
    fill <- 0
    # rowdims <- data[col__ == 1, (row.vars), with = FALSE]
    # coldims <- data[row__ == 1, (col.vars), with = FALSE]
  } else {
    # rowdims <- unique(data[, (row.vars), with = FALSE])
    # coldims <- unique(data[, (col.vars), with = FALSE])
  }
  rowdims <- unique(data[, (row.vars), with = FALSE])
  coldims <- unique(data[, (col.vars), with = FALSE])
  data.m <- matrix(fill[1], nrow = max(data[["row__"]]),
                   ncol = max(data[["col__"]]))
  data.m[cbind(data[["row__"]], data[["col__"]])] <- data[[value.var]]
  
  return(list(matrix = data.m,
              coldims = coldims,
              rowdims = rowdims))
}

makematrix <- function(data1,f){
  if (length(f) == 1) {
    f <- stringr::str_split(f, "~", n = 2)[[1]]
  }else{
    f <- f[-1]
  }
  # Define row.vars and col.vars from formula: aim: row.vars=c("lat","lon"), col.vars="date"
  value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
  matrix.vars <- f[stringr::str_detect(f, "\\|")]
  matrix.vars <- stringr::str_split(matrix.vars, "\\|", n = 2)[[1]]
  row.vars <- stringr::str_squish(stringr::str_split(matrix.vars[1],"\\+")[[1]])
  col.vars <- stringr::str_squish(stringr::str_split(matrix.vars[2],"\\+")[[1]])
  
  dcast.formula <- stringr::str_squish(f[stringr::str_detect(f,"\\|")])
  dcast.formula <- stats::as.formula(stringr::str_replace(dcast.formula,"\\|", "~"))
  value.var <- stringr::str_squish(f[!stringr::str_detect(f,"\\|")])
  
  A <- .tidy2matrix(data1, dcast.formula, value.var, fill = NULL) # 694 x 852
  
  return(A)
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
subset = list(lat = list(-46:0), lon = list(0:20,310:360))

#iyr: 1960:2017, mmb: 1:15
iyr=1960

# Only for the first member
mmb=1
# Form input file name
ifile=paste0(ipath,"tos_Omon_EC-Earth3_dcppA-hindcast_s",iyr,"-r",mmb,"i2p1f1_gn_",iyr,"11-",(iyr+10),"12_r180x91.nc")
# Open file as data.frame
data1 = metR::ReadNetCDF(ifile,vars=c("tos"), out='data.frame',subset = subset)
# Modify name of variables
colnames(data1) = c("targetdate","lat","lon",paste0("sst.",mmb))

# For the rest of the members, merge sst with the previous data.frame
for(mmb in 2:15){
  ifile=paste0(ipath,"tos_Omon_EC-Earth3_dcppA-hindcast_s",iyr,"-r",mmb,"i2p1f1_gn_",iyr,"11-",(iyr+10),"12_r180x91.nc")
  data = metR::ReadNetCDF(ifile,vars=c("tos"), out='data.frame',subset = subset)
  colnames(data) = c("targetdate","lat","lon",paste0("sst.",mmb))
  data1=merge(data1,data,by=c("targetdate","lat","lon"),all=TRUE)
}
rm(data)
as.data.table(data1)
# Add start date as a column
data1[,startdate := as.Date(as.character(paste0(iyr,"-11-01")))]

sst.ECE=data1

# Continue with other years (CHANGE END TO 2017)
for(iyr in 1961:2017){
  # Only for the first member
  mmb=1
  # Form input file name
  ifile=paste0(ipath,"tos_Omon_EC-Earth3_dcppA-hindcast_s",iyr,"-r",mmb,"i2p1f1_gn_",iyr,"11-",(iyr+10),"12_r180x91.nc")
  # Open file as data.frame
  data1 = metR::ReadNetCDF(ifile,vars=c("tos"), out='data.frame',subset = subset)
  # Modify name of variables
  colnames(data1) = c("targetdate","lat","lon",paste0("sst.",mmb))
  
  # For the rest of the members, merge tos with the previous data.frame
  for(mmb in 2:15){
    ifile=paste0(ipath,"tos_Omon_EC-Earth3_dcppA-hindcast_s",iyr,"-r",mmb,"i2p1f1_gn_",iyr,"11-",(iyr+10),"12_r180x91.nc")
    data = metR::ReadNetCDF(ifile,vars=c("tos"), out='data.frame',subset = subset)
    colnames(data) = c("targetdate","lat","lon",paste0("sst.",mmb))
    data1=merge(data1,data,by=c("targetdate","lat","lon"),all=TRUE)
  }
  rm(data)
  as.data.table(data1)
  # Add start date as a column
  data1[,startdate := as.Date(as.character(paste0(iyr,"-11-01")))]
  
  sst.ECE=rbind(sst.ECE,data1)
  rm(data1)
}

sst.ECE$targetdate=as.Date(sst.ECE$targetdate)

#________________________________________
# Computes Ensemble Mean and Climatology \______________________________________________

# Computes ensemble mean
sst.ECE[, sst.em := mean(as.matrix(.SD), na.rm = TRUE), by = .(lon, lat, targetdate,startdate)]

# # Remove linear trend
# # Order according to targetdate
# sst.ECE=sst.ECE[order(targetdate),]
# sst.ECE=sst.ECE[,sst.em := detrend(sst.em),by=.(lat,lon)]
# for (mmb in 1:15) {
#   sst.ECE[ , eval(parse(text = paste0("sst.",mmb, ":=","detrend(sst.",mmb,")" ))),by=.(lat,lon)]
# }

# Adds lead (in months, lead=1 for November)
sst.ECE[, lead := (month(targetdate)-month(startdate)+1+12*(year(targetdate)-year(startdate)))]
sst.ECE[, targetmonth := (month(targetdate))]

# Compute lead-dependent climatology
clim.sst.ECE = sst.ECE[, .(clim = ave(sst.em)), by=.(lat,lon,targetmonth,lead)]
# Remove repeated rows
clim.sst.ECE = unique(clim.sst.ECE, by=c("lat","lon","targetmonth","lead"))
# Prepare variable
clim.sst.ECE$smthclim=NA_real_
# Compute 3 point running mean
clim.sst.ECE = clim.sst.ECE[,.(targetmonth,lead,clim,smthclim = fun3mrm(clim)),by=.(lat,lon)]
clim.sst.ECE$clim=NULL

# Define lead ime in years
sst.ECE[, lead.year := ((year(targetdate)-year(startdate)))]
# Adjust december to be part of lead +1 year (year of JF)
sst.ECE[targetmonth==12,]$lead.year=sst.ECE[targetmonth==12,]$lead.year+1

#________________________________________
# Compute anomalies                      \______________________________________________

# Merge climatology
sst.ECE=merge(sst.ECE,clim.sst.ECE,by=c("lat","lon","targetmonth","lead"),all.x=TRUE)
# Move columns to the left
setcolorder(sst.ECE,neworder=c("lat","lon","startdate","targetdate","targetmonth","lead","lead.year"))

# Compute anomalies for each member
for (mmb in 1:15) {
  sst.ECE[ , eval(parse(text = paste0("asst.",mmb, ":=","sst.",mmb,"-smthclim" )))]
}
# Compute anomaly for the ensemble mean
sst.ECE[,asst.em:=sst.em-smthclim]

# Remove raw values from data table
for (mmb in 1:15) {
  eval(parse(text = paste0("sst.ECE$sst.",mmb, "=","NULL" )))
}
sst.ECE$sst.em=NULL

sst.ECE$smthclim=NULL

#________________________________________
# Save data table                        \______________________________________________

save("sst.ECE",file="/home/maralv/data/asst.ECEarth.DCP.19612017.RData")
