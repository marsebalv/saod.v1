# This script computes the SVD of sst&slp in the historical EC-Earth3 runs
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
# Loads data                              \_____________________________________________

load(file="/home/maralv/data/asst.ECEarth.hist.19492014.RData")
load(file="/home/maralv/data/apsl.ECEarth.hist.19492014.RData")

#################### Settings ########################
# Remove trend? TRUE or FALSE
remove.trend=TRUE

# Members being used
members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)
# Period selection
start.date=as.Date("1979-12-01")
end.date=as.Date("2014-12-16")
######################################################
# Select range of years
st.yr=year(start.date)
ed.yr=year(end.date)

sst.ECE=sst.ECE[targetdate>=start.date & targetdate<end.date,]
psl.ECE=psl.ECE[targetdate>=start.date & targetdate<end.date,]

# Rearrange data using number of member as variable
# STT
# rename according to member
old = c("asst.1","asst.2","asst.4","asst.5","asst.6","asst.7","asst.8","asst.9","asst.10","asst.11","asst.12","asst.13","asst.14","asst.15","asst.16","asst.em")
new = c("1","2","4","5","6","7","8","9","10","11","12","13","14","15","16","99")
setnames(sst.ECE,old,new)
# Melt data table for easy operation
sst.ECE=melt(sst.ECE,id=c("lat","lon","targetdate","targetmonth"),variable.name="member",value.name="asst")

sst.ECE$member=as.numeric(as.character(sst.ECE$member))

# SLP
# rename according to member
old = c("apsl.1","apsl.2","apsl.4","apsl.5","apsl.6","apsl.7","apsl.8","apsl.9","apsl.10","apsl.11","apsl.12","apsl.13","apsl.14","apsl.15","apsl.16","apsl.em")
setnames(psl.ECE,old,new)
# Melt data table for easy operation
psl.ECE=melt(psl.ECE,id=c("lat","lon","targetdate","targetmonth"),variable.name="member",value.name="apsl")  

psl.ECE$member=as.numeric(as.character(psl.ECE$member))

# Retain only the ensemble mean
sst.ECE=sst.ECE[member==99,]
psl.ECE=psl.ECE[member==99,]

sst.ECE$member=NULL
psl.ECE$member=NULL

#________________________________________
# Linear trend removal while monthly     \______________________________________________

if(remove.trend==TRUE){
  # Remove linear trend of SST anomalies
  sst.ECE=sst.ECE[order(targetdate),]
  sst.ECE=sst.ECE[,dt.asst := detrend(asst),by=.(lat,lon)]
  sst.ECE$asst=NULL
  setnames(sst.ECE,"dt.asst","asst")
  
  # Remove linear trend of PSL anomalies
  psl.ECE=psl.ECE[order(targetdate),]
  psl.ECE=psl.ECE[,dt.apsl := detrend(apsl),by=.(lat,lon)]
  psl.ECE$apsl=NULL
  setnames(psl.ECE,"dt.apsl","apsl")
}

#________________________________________
# Apply latitude weight to the anomalies \______________________________________________

sst.ECE$asst=sst.ECE$asst*sqrt(cos(sst.ECE$lat*pi/180))
psl.ECE$apsl=psl.ECE$apsl*sqrt(cos(psl.ECE$lat*pi/180))

#________________________________________
# Remove NA values                       \______________________________________________

sst.ECE = sst.ECE[!is.na(asst),]
psl.ECE = psl.ECE[!is.na(apsl),]

##############################################
#
#   SVD - Atlantic SST & SLP
#
##############################################

# Change data to rows=space, columns=time

# Rename targetdate
setnames(sst.ECE,"targetdate","date")
setnames(psl.ECE,"targetdate","date")

S = makematrix(sst.ECE,"asst ~ lat + lon | date")
P = makematrix(psl.ECE,"apsl ~ lat + lon | date")
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
exp.coef = data.table(date = unique(sst.ECE$date),ec.sst.1=A[,1],ec.sst.2=A[,2],ec.sst.3=A[,3],ec.slp.1=B[,1],ec.slp.2=B[,2],ec.slp.3=B[,3])

# Create data table with smoothed expansion coefficients and normalize respect to sd
exp.coef.norm = data.table(date = unique(sst.ECE$date),ec.sst.1=rollmean(A[,1],13,align="center",fill=NA),ec.sst.2=rollmean(A[,2],13,align="center",fill=NA),ec.sst.3=rollmean(A[,3],13,align="center",fill=NA),ec.slp.1=rollmean(B[,1],13,align="center",fill=NA),ec.slp.2=rollmean(B[,2],13,align="center",fill=NA),ec.slp.3=rollmean(B[,3],13,align="center",fill=NA))
#exp.coef.norm = data.table(date = unique(sst.ECE$date),ec.sst.1=A[,1],ec.sst.2=A[,2],ec.sst.3=A[,3],ec.slp.1=B[,1],ec.slp.2=B[,2],ec.slp.3=B[,3])
exp.coef.norm = exp.coef.norm[,.(date,ec.sst.1=ec.sst.1/sd(ec.sst.1,na.rm=TRUE),ec.sst.2=ec.sst.2/sd(ec.sst.2,na.rm=TRUE),ec.sst.3=ec.sst.3/sd(ec.sst.3,na.rm=TRUE),ec.slp.1=ec.slp.1/sd(ec.slp.1,na.rm=TRUE),ec.slp.2=ec.slp.2/sd(ec.slp.2,na.rm=TRUE),ec.slp.3=ec.slp.3/sd(ec.slp.3,na.rm=TRUE))]

rm(A,B,C,SV,SV2,S,P)

# As I will present spatial patterns as homogeneous correlation maps:
psl.ECE = merge(psl.ECE,exp.coef[,.(date,ec.slp.1,ec.slp.2,ec.slp.3)],by="date")
sst.ECE = merge(sst.ECE,exp.coef[,.(date,ec.sst.1,ec.sst.2,ec.sst.3)],by="date")

psl.ECE = psl.ECE[,.(date,apsl,ec.slp.1,ec.slp.2,ec.slp.3,hcm.1 = cor(apsl,ec.slp.1),hcm.2 = cor(apsl,ec.slp.2),hcm.3 = cor(apsl,ec.slp.3)),by=.(lat,lon)]
sst.ECE = sst.ECE[,.(date,asst,ec.sst.1,ec.sst.2,ec.sst.3,hcm.1 = cor(asst,ec.sst.1),hcm.2 = cor(asst,ec.sst.2),hcm.3 = cor(asst,ec.sst.3)),by=.(lat,lon)]

#________________________________________
# Plotting and saving                    \______________________________________________

# Change longitude to -180:180 to plot correctly
psl.ECE[lon>180]$lon=psl.ECE[lon>180]$lon-360
sst.ECE[lon>180]$lon=sst.ECE[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-0.9
bmax=0.9
bbreaks.contours=c(-99,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,99)

# Need to select a single date as values are repeated by date
g1 = ggplot() +
  geom_contour_fill(data=sst.ECE[date==sst.ECE$date[1]],aes(lon, lat, z = hcm.1, fill=stat(level)),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       limits=c(bmin,bmax),
                       super = ScaleDiscretised,
                       guide = guide_colorsteps(),
                       oob  = scales::squish)+
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 10))+
  new_scale_color() +
  geom_contour(data=psl.ECE[date==psl.ECE$date[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=psl.ECE[date==psl.ECE$date[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("SVD1 as homogeneous correlation map: SST (shaded), SLP (contours) SCF: ",as.character(round(SCF[1],1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.norm,aes(date, ec.sst.1),col="red")+
  geom_line(data=exp.coef.norm,aes(date, ec.slp.1),col="black")+
  scale_x_date(breaks=exp.coef.norm$date[seq(1,length(exp.coef.norm$date),48)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-5,5,1),limits=c(-6,6),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("Normalized Expansion Coefficients: SST (red), SLP (black). r=",as.character(round(cor(exp.coef$ec.sst.1,exp.coef$ec.slp.1,use = "pairwise.complete.obs"),2))," for raw EC"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

# Save Figure and PCs for the lead

if(remove.trend==TRUE){

    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0("All months , historical: SVD of SST-SLP anomalies (no trend, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/AllMonths_historical_SVD_SST_SLP_weighted_notrend_ensmean_",st.yr,"-",ed.yr,".png"),plot=fig,width = 8, height = 8)
    
  
}else{
  # save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/AllMonths_historical_ECEarth3_PCs_EOFs_SST_weighted_allmembers.RData"))

    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0("All months, historical: SVD of SST-SLP anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/AllMonths_historical_SVD_SST_SLP_weighted_ensmean_",st.yr,"-",ed.yr,".png"),plot=fig,width = 8, height = 8)

} 

