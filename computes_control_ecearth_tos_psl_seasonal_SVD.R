# This script computes the SVD of sst&slp in the control EC-Earth3 run
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

load(file="/home/maralv/data/asst.ECEarth.control.18502849.RData")
load(file="/home/maralv/data/apsl.ECEarth.control.18502849.RData")

#################### Settings ########################
# Remove trend? TRUE or FALSE
remove.trend=FALSE
# Select season
sel.season="DJF"
# Rolling years
roll.years=FALSE
# number of years to roll
ny=4
# Period selection
start.date=as.Date("1850-12-01")
end.date=as.Date("2849-12-16")
######################################################
# Select range of years
st.yr=year(start.date)
ed.yr=year(end.date)

sst.ECE=sst.ECE[date>=start.date & date<end.date,]
psl.ECE=psl.ECE[date>=start.date & date<end.date,]

#________________________________________
# Linear trend removal while monthly     \______________________________________________

if(remove.trend==TRUE){
  # Remove linear trend of SST anomalies
  sst.ECE=sst.ECE[order(date),]
  sst.ECE=sst.ECE[,dt.asst := detrend(asst),by=.(lat,lon)]
  sst.ECE$asst=NULL
  setnames(sst.ECE,"dt.asst","asst")
  
  # Remove linear trend of PSL anomalies
  psl.ECE=psl.ECE[order(date),]
  psl.ECE=psl.ECE[,dt.apsl := detrend(apsl),by=.(lat,lon)]
  psl.ECE$apsl=NULL
  setnames(psl.ECE,"dt.apsl","apsl")
}

#________________________________________
# Perform seasonal averages              \______________________________________________

# SST

# Define seasons
sst.ECE = sst.ECE[month==12 | month==1 | month==2, season := "DJF" ]
sst.ECE = sst.ECE[month==3 | month==4 | month==5, season := "MAM" ]
sst.ECE = sst.ECE[month==6 | month==7 | month==8, season := "JJA" ]
sst.ECE = sst.ECE[month==9 | month==10 | month==11, season := "SON" ]

sst.ECE = sst.ECE[season == sel.season,]

# Add year for season
sst.ECE$s.year=year(sst.ECE$date)
sst.ECE[month==12,]$s.year=sst.ECE[month==12,]$s.year+1

# Perform seasonal averages
sst.ECE=sst.ECE[,s.asst := ave(asst),by=.(s.year,lat,lon)]

# Eliminate monthly data & clean data table
sst.ECE$asst=NULL
sst.ECE$month=NULL
sst.ECE$season=NULL

# Remove repeated rows
sst.ECE = unique(sst.ECE, by=c("lat","lon","s.year"))

# Compute seasonal average for the ensemble mean
s.clim.sst.ECE=sst.ECE[,.(lat,lon,s.year,s.asst)]
s.clim.sst.ECE=s.clim.sst.ECE[,s.clim := ave(s.asst),by=.(lat,lon)]
s.clim.sst.ECE=s.clim.sst.ECE[,.(lat,lon,s.clim)]
s.clim.sst.ECE=unique(s.clim.sst.ECE)
# Merge and compute seasonal anomaly
sst.ECE = merge(sst.ECE,s.clim.sst.ECE)
# Compute anomaly
sst.ECE=sst.ECE[,sa.asst := (s.asst - s.clim)]
sst.ECE$s.asst=NULL
sst.ECE$s.clim=NULL
setnames(sst.ECE,"sa.asst","s.asst")

# SLP
# Define seasons
psl.ECE = psl.ECE[month==12 | month==1 | month==2, season := "DJF" ]
psl.ECE = psl.ECE[month==3 | month==4 | month==5, season := "MAM" ]
psl.ECE = psl.ECE[month==6 | month==7 | month==8, season := "JJA" ]
psl.ECE = psl.ECE[month==9 | month==10 | month==11, season := "SON" ]

psl.ECE = psl.ECE[season == sel.season,]

# Add year for season
psl.ECE$s.year=year(psl.ECE$date)
psl.ECE[month==12,]$s.year=psl.ECE[month==12,]$s.year+1


# Perform seasonal averages
psl.ECE=psl.ECE[,s.apsl := ave(apsl),by=.(s.year,lat,lon)]

# Eliminate monthly data & clean data table
psl.ECE$apsl=NULL
psl.ECE$month=NULL
psl.ECE$season=NULL

# Remove repeated rows
psl.ECE = unique(psl.ECE, by=c("lat","lon","s.year"))

# Compute seasonal average for the ensemble mean
s.clim.psl.ECE=psl.ECE[,.(lat,lon,s.year,s.apsl)]
s.clim.psl.ECE=s.clim.psl.ECE[,s.clim := ave(s.apsl),by=.(lat,lon)]
s.clim.psl.ECE=s.clim.psl.ECE[,.(lat,lon,s.clim)]
s.clim.psl.ECE=unique(s.clim.psl.ECE)
# Merge and compute seasonal anomaly
psl.ECE = merge(psl.ECE,s.clim.psl.ECE)
# Compute anomaly
psl.ECE=psl.ECE[,sa.apsl := (s.apsl - s.clim)]
psl.ECE$s.apsl=NULL
psl.ECE$s.clim=NULL
setnames(psl.ECE,"sa.apsl","s.apsl")

#________________________________________
# Perform ny-year rolling means           \______________________________________________

if(roll.years==TRUE){
  # SST
  sst.ECE=sst.ECE[order(date),]
  alig="left" #to be used for the rolling mean, using for time=t data between t and t+3
  
  sst.ECE=sst.ECE[, s.year.roll := frollmean(s.year,n=ny,align=alig),by=.(lat,lon)]
  sst.ECE=sst.ECE[, s.asst.roll := frollmean(s.asst,n=ny,align=alig),by=.(lat,lon)]
  
  # Remove rows which could not be used to compute the rolling mean
  sst.ECE = sst.ECE[!is.na(s.year.roll),]
  
  # Rearrange variable names
  sst.ECE$s.asst=NULL
  setnames(sst.ECE,"s.asst.roll","asst")
  
  #PSL
  psl.ECE=psl.ECE[order(date),]
  alig="left" #to be used for the rolling mean, using for time=t data between t and t+3
  
  psl.ECE=psl.ECE[, s.year.roll := frollmean(s.year,n=ny,align=alig),by=.(lat,lon)]
  psl.ECE=psl.ECE[, s.apsl.roll := frollmean(s.apsl,n=ny,align=alig),by=.(lat,lon)]
  
  # Remove rows which could not be used to compute the rolling mean
  psl.ECE = psl.ECE[!is.na(s.year.roll),]
  
  # Rearrange variable names
  psl.ECE$s.apsl=NULL
  setnames(psl.ECE,"s.apsl.roll","apsl")
}else{
  setnames(sst.ECE,"s.asst","asst")
  setnames(psl.ECE,"s.apsl","apsl")
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
# exp.coef.norm = data.table(date = unique(sst.ECE$date),ec.sst.1=rollmean(A[,1],13,align="center",fill=NA),ec.sst.2=rollmean(A[,2],13,align="center",fill=NA),ec.sst.3=rollmean(A[,3],13,align="center",fill=NA),ec.slp.1=rollmean(B[,1],13,align="center",fill=NA),ec.slp.2=rollmean(B[,2],13,align="center",fill=NA),ec.slp.3=rollmean(B[,3],13,align="center",fill=NA))
exp.coef.norm = data.table(date = unique(sst.ECE$date),ec.sst.1=A[,1],ec.sst.2=A[,2],ec.sst.3=A[,3],ec.slp.1=B[,1],ec.slp.2=B[,2],ec.slp.3=B[,3])
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

bmin=-0.8
bmax=0.8
bbreaks.contours=c(-99,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,99)

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
  scale_x_date(breaks=exp.coef.norm$date[seq(1,length(exp.coef.norm$date),30)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-5,5,1),limits=c(-6,6),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("Normalized Expansion Coefficients: SST (red), SLP (black). r=",as.character(round(cor(exp.coef$ec.sst.1,exp.coef$ec.slp.1,use = "pairwise.complete.obs"),2))," for raw EC"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

# Save Figure and PCs for the lead

if(remove.trend==TRUE){
  if(roll.years==TRUE){
    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0(sel.season," , control: SVD of SST-SLP anomalies (no trend, rolling ",ny," years, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_control_SVD_SST_SLP_weighted_notrend_rollingyears_",ny,"_",st.yr,"-",ed.yr,".png"),plot=fig,width = 8, height = 8)
    
  }else{
    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0(sel.season," , control: SVD of SST-SLP anomalies (no trend, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_control_SVD_SST_SLP_weighted_notrend_",st.yr,"-",ed.yr,".png"),plot=fig,width = 8, height = 8)
    
  }
  
}else{
  # save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/",sel.season,"_control_ECEarth3_PCs_EOFs_SST_weighted_allmembers.RData"))
  if(roll.years==TRUE){
    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0(sel.season," , control: SVD of SST-SLP anomalies (rolling ",ny," years, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_control_SVD_SST_SLP_weighted_rollingyears_",ny,"_",st.yr,"-",ed.yr,".png"),plot=fig,width = 8, height = 8)
    
  }else{
    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0(sel.season," , control: SVD of SST-SLP anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_control_SVD_SST_SLP_weighted_",st.yr,"-",ed.yr,".png"),plot=fig,width = 8, height = 8)
    
  }
} 

