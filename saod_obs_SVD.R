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
#Datos de SLP (NCEP-NCAR Reana1)
# lat e/2.5°
# lon e/2.5°

# File to load
file.in="/home/maralv/data/slp.ncepncar.mon.nc"

# Longitude regions (divided in the Atlantic by Greenwich)
lon1=seq(0,20,2.5)
lon2=seq(310,357.5,2.5)

subset.in1 = list(Y = seq(-45,0,2.5), X=lon1)
subset.in2 = list(Y = seq(-45,0,2.5), X=lon2)
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

# Ojo, no tengo que pasarle como argumento a la funci+on "clim$clim1", y tengo que repetir,
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

# Change longitude to -180:180 to plot correctly
sst[lon>180]$lon=sst[lon>180]$lon-360
slp[lon>180]$lon=slp[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

# Weighted by latitude
sst$ssta=sst$ssta*sqrt(cos(sst$lat*pi/180))
slp$slpa=slp$slpa*sqrt(cos(slp$lat*pi/180))

##############################################
#
#   SVD - Atlantic SST & SLP
#
##############################################

# Change data to rows=space, columns=time

S = makematrix(sst,"ssta ~ lat + lon | date")
P = makematrix(slp,"slpa ~ lat + lon | date")

# Perform SVD

nmodes=3
nt=length(P$coldims$date)
# S: locations x times
# P: locations x times
C = (S$matrix %*% t(P$matrix))/(nt - 1) # covariance matrix
# Perform SVD
MCA = svd(C) #C = U D V'
# Columns of U are singular vectors of S 
# Columns of V are singular vectors of P

# Expansion Coefficients
A = t(S$matrix) %*% MCA$u #(la tenia arreglada distinto, es times x loc.sst %*% loc.sst x loc.slp )
B = t(P$matrix) %*% MCA$v

# Squared Covariance and Squared Covariance Fraction
SV = MCA$d
SV2 = SV * SV
S = sum(SV2)
SCF = SV2/S
nm = min(nmodes, length(SV))
SCF = SCF[1:nm] * 100

# Create data table with expansion coefficients
exp.coef = data.table(date = unique(sst$date),ec.sst.1=A[,1],ec.sst.2=A[,2],ec.sst.3=A[,3],ec.slp.1=B[,1],ec.slp.2=B[,2],ec.slp.3=B[,3])

# Create data table with smoothed expansion coefficients and normalize respect to sd
exp.coef.norm = data.table(date = unique(sst$date),ec.sst.1=rollmean(A[,1],13,align="center",fill=NA),ec.sst.2=rollmean(A[,2],13,align="center",fill=NA),ec.sst.3=rollmean(A[,3],13,align="center",fill=NA),ec.slp.1=rollmean(B[,1],13,align="center",fill=NA),ec.slp.2=rollmean(B[,2],13,align="center",fill=NA),ec.slp.3=rollmean(B[,3],13,align="center",fill=NA))
exp.coef.norm = exp.coef.norm[,.(date,ec.sst.1=ec.sst.1/sd(ec.sst.1,na.rm=TRUE),ec.sst.2=ec.sst.2/sd(ec.sst.2,na.rm=TRUE),ec.sst.3=ec.sst.3/sd(ec.sst.3,na.rm=TRUE),ec.slp.1=ec.slp.1/sd(ec.slp.1,na.rm=TRUE),ec.slp.2=ec.slp.2/sd(ec.slp.2,na.rm=TRUE),ec.slp.3=ec.slp.3/sd(ec.slp.3,na.rm=TRUE))]

rm(A,B,C,SV,SV2,S,P)

# As I will present spatial patterns as homogeneous correlation maps:
slp = merge(slp,exp.coef[,.(date,ec.slp.1,ec.slp.2,ec.slp.3)],by="date")
sst = merge(sst,exp.coef[,.(date,ec.sst.1,ec.sst.2,ec.sst.3)],by="date")

slp = slp[,.(date,slpa,ec.slp.1,ec.slp.2,ec.slp.3,hcm.1 = cor(slpa,ec.slp.1),hcm.2 = cor(slpa,ec.slp.2),hcm.3 = cor(slpa,ec.slp.3)),by=.(lat,lon)]
sst = sst[,.(date,ssta,ec.sst.1,ec.sst.2,ec.sst.3,hcm.1 = cor(ssta,ec.sst.1),hcm.2 = cor(ssta,ec.sst.2),hcm.3 = cor(ssta,ec.sst.3)),by=.(lat,lon)]

# Change longitude to -180:180 to plot correctly
sst[lon>180]$lon=sst[lon>180]$lon-360
slp[lon>180]$lon=slp[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

# Only one plot, both patterns
bmin=-0.9
bmax=0.9
bbreaks=c(-99,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,99)

# Need to select a single date as values are repeated by date
g1 = ggplot() +
  geom_contour_fill(data=sst[date==as.Date("1949-01-01")],aes(lon, lat, z = hcm.1, fill=stat(level)),breaks=bbreaks,na.fill=TRUE)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       limits=c(bmin,bmax),
                       super = ScaleDiscretised,
                       guide = guide_colorsteps(),
                       oob  = scales::squish)+
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 10))+
  new_scale_color() +
  geom_contour(data=slp[date==as.Date("1949-01-01")],aes(lon, lat, z = hcm.1),bbreaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=slp[date==as.Date("1949-01-01")],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("SVD1 as homogeneous correlation map: SST (shaded), SLP (contours) SCF: ",as.character(round(SCF[1],1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.norm,aes(date, ec.sst.1),col="red")+
  geom_line(data=exp.coef.norm,aes(date, ec.slp.1),col="black")+
  scale_x_date(breaks=exp.coef.norm$date[seq(1,length(exp.coef.norm$date),48)],date_labels = "%Y-%m",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.5,3.5),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("Smoothed Expansion Coefficients: SST (red), SLP (black). r=",as.character(round(cor(exp.coef$ec.sst.1,exp.coef$ec.slp.1,use = "pairwise.complete.obs"),2))," for raw EC"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0("Leading SVD (SLP,SST) (without trend, weighted by cos(lat)) following Venegas et al. (1997)"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/Leading_SVD_SSTA_SLPA_weighted_notrend.png"),plot=fig,width = 8, height = 8)

save("exp.coef","exp.coef.norm",file="SAO_monthly_expansion_coefs.RData")
