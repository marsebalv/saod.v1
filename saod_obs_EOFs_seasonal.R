# This script opens data and performs SVD to obtain first modes of South Atlantic
# MSLP and SST.
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

library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#---------------------------------------------------------------------------------------
#  Functions to compute power spectra by Elio Campitelli (GitHub @eliocamp)
#---------------------------------------------------------------------------------------
#' FFT spectrum
#'
#' A light wrapper around [stats::spec.pgram()] and [stats::ar()]. Returns the spectrum of
#' the signal and the spectrum of the fitted autoregressiv model.
#'
#' @param x numeric vector
#' @param B number of bootstrap samples.
#' @param spans vector of odd integers giving the widths of modified Daniell
#' smoothers to be used to smooth the periodogram
#' @param ... other arguments passed to [stats::spec.pgram()]
#'
#' @export
fftspectrum <- function(x, spans = NULL, B = 10000, ..., probs = 0.95) {
  mtm <- spec.pgram(x, spans = spans, ..., plot = FALSE)  
  out <- as.data.table(mtm[c("freq", "spec")])  
  ar <- ar(ts(x))
  # rho <- a$ar
  # var <- a$var.pred
  out[, ar_spectrum := arspectrum(mtm$freq, ar$ar, ar$var.pred)]
  out[, c(scales::percent(probs)) := null_ar_spectrum(B = B, length(x), ar, spans = spans, ..., probs = probs)]  
  return(out[])
}

arspectrum <- function(freq, rho, var) {
  k <- seq_along(rho)
  e <- vapply(freq, function(f) sum(rho * exp(-2*1i*pi*f*k)), complex(1))
  var / (Mod(1 - e))^2
}

null_ar_spectrum_ <- function(B = 100, n, ar, spans = NULL, ..., probs = 0.95) {
  y <- as.vector(arima.sim(model = list(ar = ar$ar), n = n))*sqrt(ar$var.pred)
  nfreq <- length(spec.pgram(y, spans = spans, ..., plot = FALSE)$spec)
  boots <- vapply(seq_len(B), function(b) {
    y <- as.vector(arima.sim(model = list(ar = ar$ar), n = n))*sqrt(ar$var.pred)
    spec.pgram(y, spans = spans, plot = FALSE)$spec  
    
  }, numeric(nfreq))
  data.table::as.data.table((apply(boots, 1, quantile, probs = probs)))
}

null_ar_spectrum <- memoise::memoise(null_ar_spectrum_)
# null_ar_spectrum <- null_ar_spectrum_

null_spec <- memoise::memoise(function(x, spans, B = 1000, ..., probs = 0.95) {
  b <- boot::boot(x, function(d, i)  spec.pgram(d[i], spans = spans,
                                                ...,
                                                plot = FALSE)$spec,
                  R = B)  
  apply(b$t, 2, quantile, probs = probs)
})

#---------------------------------------------------------------------------------------
#  Man Program
#---------------------------------------------------------------------------------------

#################### Settings ########################
# Remove trend? TRUE or FALSE
remove.trend=TRUE
# Select season
sel.season="DJF"
######################################################

if(sel.season == "DJF"){
  mmdd="-01-01"
}else if(sel.season == "MAM"){
  mmdd="-03-01"
}else if(sel.season == "JJA"){
  mmdd="-06-01"
}else if(sel.season == "SON"){
  mmdd="-09-01"
}

#Datos de SLP (NCEP-NCAR Reana1)
# lat e/2.5°
# lon e/2.5°

# File to load
file.in="/home/maralv/data/slp.ncepncar.mon.nc"

# Longitude regions (divided in the Atlantic by Greenwich)
lon1=seq(0,20,2.5)
lon2=seq(310,357.5,2.5)

subset.in1 = list(Y = seq(-46,0,2.5), X=lon1)
subset.in2 = list(Y = seq(-46,0,2.5), X=lon2)
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
lon2=seq(310,358,2)

subset.in1 = list(Y = seq(-46,0,2), X=lon1)
subset.in2 = list(Y = seq(-46,0,2), X=lon2)

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

# No tengo que pasarle como argumento a la función "clim$clim1", y tengo que repetir,
# más allá del "by", qué otras variables quiero seguir teniendo en la DT

clim.slp = clim.slp[,.(month,smthclim = fun3mrm(clim)),by=.(lat,lon)]
clim.sst = clim.sst[,.(month,smthclim = fun3mrm(clim)),by=.(lat,lon)]

rm(mask,mask2,subset.in1,subset.in2,file.in,lon1,lon2)

# Merge climatology with values to compute anomalies
slp = merge(slp,clim.slp,by=c("lat","lon","month"),all=TRUE)
sst = merge(sst,clim.sst,by=c("lat","lon","month"),all=TRUE)

rm(clim.slp,clim.sst)

slp$slpa=slp$slp-slp$smthclim
sst$ssta=sst$sst-sst$smthclim

sst=sst[,.(lat,lon,date,ssta)]
slp=slp[,.(lat,lon,date,slpa)]

# sst has all missing values only for lon,lat=(322,-4) and all dates. Removing the data point.
sst=sst[!is.na(ssta),]

if(remove.trend==TRUE){
  # Remove linear trend of SST anomalies
  sst=sst[order(date),]
  sst=sst[,dtssta := detrend(ssta),by=.(lat,lon)]
  sst$ssta=NULL
  setnames(sst,"dtssta","ssta")

  # Remove linear trend of SLP anomalies
  slp=slp[order(date),]
  slp=slp[,dtslpa := detrend(slpa),by=.(lat,lon)]
  slp$slpa=NULL
  setnames(slp,"dtslpa","slpa")
}


############
# SST season selection
###########
# Define seasons
sst = sst[month(sst$date)==12 | month(sst$date)==1 | month(sst$date)==2, season := "DJF" ]
sst = sst[month(sst$date)==3 | month(sst$date)==4 | month(sst$date)==5, season := "MAM" ]
sst = sst[month(sst$date)==6 | month(sst$date)==7 | month(sst$date)==8, season := "JJA" ]
sst = sst[month(sst$date)==9 | month(sst$date)==10 | month(sst$date)==11, season := "SON" ]

sst = sst[,season.year := year(date)]
# Adjust december to be part of season.year +1 year (year of JF)
sst[month(date)==12,]$season.year=sst[month(date)==12,]$season.year+1

sst = sst[season == sel.season,]

# Perform seasonal averages
sst = sst[,s.ssta := ave(ssta),by=.(season.year,lat,lon)]

# Eliminate monthly data
sst$ssta=NULL

# Clean DT
sst$season=NULL
sst$date=NULL

# Remove repeated rows
sst = unique(sst, by=c("lat","lon","season.year"))

# rename
setnames(sst,"s.ssta","ssta")
setnames(sst,"season.year","date")

# Modify date just to configure variable as date (only year is valid)
sst[,date := as.Date(as.character(paste0(date,mmdd)))]

############
# SLP season selection
###########

# Define seasons
slp = slp[month(slp$date)==12 | month(slp$date)==1 | month(slp$date)==2, season := "DJF" ]
slp = slp[month(slp$date)==3 | month(slp$date)==4 | month(slp$date)==5, season := "MAM" ]
slp = slp[month(slp$date)==6 | month(slp$date)==7 | month(slp$date)==8, season := "JJA" ]
slp = slp[month(slp$date)==9 | month(slp$date)==10 | month(slp$date)==11, season := "SON" ]

slp = slp[,season.year := year(date)]
# Adjust december to be part of season.year +1 year (year of JF)
slp[month(date)==12,]$season.year=slp[month(date)==12,]$season.year+1

slp = slp[season == sel.season,]

# Perform seasonal averages
slp = slp[,s.slpa := ave(slpa),by=.(season.year,lat,lon)]

# Eliminate monthly data
slp$slpa=NULL

# Clean DT
slp$season=NULL
slp$date=NULL

# Remove repeated rows
slp = unique(slp, by=c("lat","lon","season.year"))

# rename
setnames(slp,"s.slpa","slpa")
setnames(slp,"season.year","date")

# Modify date just to configure variable as date (only year is valid)
slp[,date := as.Date(as.character(paste0(date,mmdd)))]




# Change longitude to -180:180 to plot correctly
sst[lon>180]$lon=sst[lon>180]$lon-360
slp[lon>180]$lon=slp[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

# To recreate the patterns as shown in Venegas et al 1997 but weighted by latitude
sst$ssta=sst$ssta*sqrt(cos(sst$lat*pi/180))
slp$slpa=slp$slpa*sqrt(cos(slp$lat*pi/180))

##############################################
#
#   SST
#
##############################################
# Computes EOFs 1 to 3 without additional latitude weight (itś already weighted in line 251)
eof.sst = metR::EOF(ssta ~ lat + lon | date, n=1:3, data = sst, B = 100, probs = c(low = 0.1, hig = 0.9))

eof.sst.1 = cut(eof.sst, 1)
eof.sst.2 = cut(eof.sst, 2)
eof.sst.3 = cut(eof.sst, 3)


# Normalize respect to standard deviation
eof.sst.1$right$normpc=eof.sst.1$right$ssta/sd(eof.sst.1$right$ssta,na.rm=TRUE)
eof.sst.2$right$normpc=eof.sst.2$right$ssta/sd(eof.sst.2$right$ssta,na.rm=TRUE)
eof.sst.3$right$normpc=eof.sst.3$right$ssta/sd(eof.sst.3$right$ssta,na.rm=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

sst=merge(sst,eof.sst.1$right[,.(date,normpc)],by=c("date"),all=TRUE)
setnames(sst,"normpc","normpc1")
sst=merge(sst,eof.sst.2$right[,.(date,normpc)],by=c("date"),all=TRUE)
setnames(sst,"normpc","normpc2")
sst=merge(sst,eof.sst.3$right[,.(date,normpc)],by=c("date"),all=TRUE)
setnames(sst,"normpc","normpc3")

sst=sst[,.(date,ssta,normpc1,normpc2,normpc3,hcorrmap.eof1=cor(ssta,normpc1, use = "pairwise.complete.obs")),by=.(lat,lon)]
sst=sst[,.(date,ssta,normpc1,normpc2,normpc3,hcorrmap.eof1,hcorrmap.eof2=cor(ssta,normpc2, use = "pairwise.complete.obs")),by=.(lat,lon)]
sst=sst[,.(date,ssta,normpc1,normpc2,normpc3,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3=cor(ssta,normpc3, use = "pairwise.complete.obs")),by=.(lat,lon)]

# Separate data tables for eofs and pcs
sst.eof=unique(sst[,.(lat,lon,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3)],by=c("lat","lon"))
sst.pcs=unique(sst[,.(date,normpc1,normpc2,normpc3)],by=c("date"))

# Plot EOFs and smoothed PCs

bmin=-0.6
bmax=0.6
bstep=0.1
bbreaks.contours=c(-99,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,99)
bbreaks.cbar=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6)
labels.cbar=as.character(bbreaks.cbar)
labels.cbar[1]=""
labels.cbar[length(labels.cbar)]=""

g1 = ggplot() +
  geom_contour_fill(data=sst.eof,aes(lon, lat, z = hcorrmap.eof1),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF1",palette="RdBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF1. Explained variance: ",as.character(round(eof.sst.1$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  geom_contour_fill(data=sst.eof,aes(lon, lat, z = hcorrmap.eof2),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF2",palette="RdBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF2. Explained variance: ",as.character(round(eof.sst.2$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g3 = ggplot() +
  geom_contour_fill(data=sst.eof,aes(lon, lat, z = hcorrmap.eof3),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF3",palette="RdBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF3. Explained variance: ",as.character(round(eof.sst.3$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g4 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(date, normpc1))+
  scale_x_date(breaks=sst.pcs$date[seq(1,length(sst.pcs$date),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC1 (normalized units)")+
  theme(axis.title = element_text(size=11))
  
g5 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(date, normpc2))+
  scale_x_date(breaks=sst.pcs$date[seq(1,length(sst.pcs$date),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC2 (normalized units)")+
  theme(axis.title = element_text(size=11))

g6 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(date, normpc3))+
  scale_x_date(breaks=sst.pcs$date[seq(1,length(sst.pcs$date),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC3 (normalized units)")+
  theme(axis.title = element_text(size=11))

if(remove.trend==TRUE){
  fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," Leading EOFs of SST anomalies (without trend, weighted by cos(lat)) following Venegas et al. (1997)"),gp=gpar(fontsize=13,font=3)))
  ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_EOFs_SSTA_weighted_notrend.png"),plot=fig,width = 10, height = 8)
}else{
  fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," Leading EOFs of SST anomalies (with trend, weighted by cos(lat)) following Venegas et al. (1997)"),gp=gpar(fontsize=13,font=3)))
  ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_EOFs_SSTA_weighted.png"),plot=fig,width = 10, height = 8)
} 


# Clean workspace
rm(fig,g1,g2,g3,g4,g5,g6)

##############################################
#
#   SLP
#
##############################################
# Computes EOFs 1 to 3 without additional latitude weight (already weighted in line 251)
eof.slp = metR::EOF(slpa ~ lat + lon | date, n=1:3, data = slp, B = 100, probs = c(low = 0.1, hig = 0.9))

eof.slp.1 = cut(eof.slp, 1)
eof.slp.2 = cut(eof.slp, 2)
eof.slp.3 = cut(eof.slp, 3)

# Normalize respect to standard deviation
eof.slp.1$right$normpc=eof.slp.1$right$slpa/sd(eof.slp.1$right$slpa,na.rm=TRUE)
eof.slp.2$right$normpc=eof.slp.2$right$slpa/sd(eof.slp.2$right$slpa,na.rm=TRUE)
eof.slp.3$right$normpc=eof.slp.3$right$slpa/sd(eof.slp.3$right$slpa,na.rm=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

slp=merge(slp,eof.slp.1$right[,.(date,normpc)],by=c("date"),all=TRUE)
setnames(slp,"normpc","normpc1")
slp=merge(slp,eof.slp.2$right[,.(date,normpc)],by=c("date"),all=TRUE)
setnames(slp,"normpc","normpc2")
slp=merge(slp,eof.slp.3$right[,.(date,normpc)],by=c("date"),all=TRUE)
setnames(slp,"normpc","normpc3")

slp=slp[,.(date,slpa,normpc1,normpc2,normpc3,hcorrmap.eof1=cor(slpa,normpc1, use = "pairwise.complete.obs")),by=.(lat,lon)]
slp=slp[,.(date,slpa,normpc1,normpc2,normpc3,hcorrmap.eof1,hcorrmap.eof2=cor(slpa,normpc2, use = "pairwise.complete.obs")),by=.(lat,lon)]
slp=slp[,.(date,slpa,normpc1,normpc2,normpc3,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3=cor(slpa,normpc3, use = "pairwise.complete.obs")),by=.(lat,lon)]

# Separate data tables for eofs and pcs
slp.eof=unique(slp[,.(lat,lon,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3)],by=c("lat","lon"))
slp.pcs=unique(slp[,.(date,normpc1,normpc2,normpc3)],by=c("date"))

# Plotting

bminbmin=-0.6
bmax=0.6
bstep=0.1
bbreaks.contours=c(-99,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,99)
bbreaks.cbar=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6)
labels.cbar=as.character(bbreaks.cbar)
labels.cbar[1]=""
labels.cbar[length(labels.cbar)]=""

g1 = ggplot() +
  geom_contour_fill(data=slp.eof,aes(lon, lat, z = hcorrmap.eof1),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF1",palette="RdYlBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF1. Explained variance: ",as.character(round(eof.slp.1$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  geom_contour_fill(data=slp.eof,aes(lon, lat, z = hcorrmap.eof2),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF2",palette="RdYlBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF2. Explained variance: ",as.character(round(eof.slp.2$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g3 = ggplot() +
  geom_contour_fill(data=slp.eof,aes(lon, lat, z = hcorrmap.eof3),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF3",palette="RdYlBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF3. Explained variance: ",as.character(round(eof.slp.3$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g4 = ggplot() +
  theme_bw()+
  geom_line(data=slp.pcs,aes(date, normpc1))+
  scale_x_date(breaks=slp.pcs$date[seq(1,length(slp.pcs$date),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC1 (normalized units)")+
  theme(axis.title = element_text(size=11))

g5 = ggplot() +
  theme_bw()+
  geom_line(data=slp.pcs,aes(date, normpc2))+
  scale_x_date(breaks=slp.pcs$date[seq(1,length(slp.pcs$date),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC2 (normalized units)")+
  theme(axis.title = element_text(size=11))

g6 = ggplot() +
  theme_bw()+
  geom_line(data=slp.pcs,aes(date, normpc3))+
  scale_x_date(breaks=slp.pcs$date[seq(1,length(slp.pcs$date),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC3 (normalized units)")+
  theme(axis.title = element_text(size=11))

if(remove.trend==TRUE){
  fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," Leading EOFs of SLP anomalies (without trend, weighted by cos(lat)) following Venegas et al. (1997)"),gp=gpar(fontsize=13,font=3)))
  ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_EOFs_SLPA_weighted_notrend.png"),plot=fig,width = 10, height = 8)
  save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/",sel.season,"_obs_PCs_EOFs_SST_weighted_notrend.RData"))
}else{
  fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," Leading EOFs of SLP anomalies (with trend, weighted by cos(lat)) following Venegas et al. (1997)"),gp=gpar(fontsize=13,font=3)))
  ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_EOFs_SLPA_weighted.png"),plot=fig,width = 10, height = 8)
  save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/",sel.season,"_obs_PCs_EOFs_SST_weighted.RData"))  
} 

# Clean workspace
rm(fig,g1,g2,g3,g4,g5,g6)


