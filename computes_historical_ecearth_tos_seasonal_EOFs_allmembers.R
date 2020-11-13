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
# Loads data                              \_____________________________________________

load(file="/home/maralv/data/asst.ECEarth.hist.19492014.RData")

#################### Settings ########################
# Remove trend? TRUE or FALSE
remove.trend=FALSE
# Select season
sel.season="DJF"
# Rolling years
roll.years=TRUE
# number of years to roll
ny=4
# Members being used
members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)
######################################################

#________________________________________
# Linear trend removal while monthly     \______________________________________________

# For now, only for ensemble mean respect to itself

if(remove.trend==TRUE){
  # Remove linear trend of SST anomalies
  sst.ECE=sst.ECE[order(targetdate),]
  sst.ECE=sst.ECE[,dt.asst.em := detrend(asst.em),by=.(lat,lon)]
  sst.ECE$asst.em=NULL
  setnames(sst.ECE,"dt.asst.em","asst.em")
  
  for (mmb in members) {
  sst.ECE[ , eval(parse(text = paste0("asst.",mmb, ":= detrend(asst.",mmb,")"))),by=.(lat,lon)]
  # eval(parse(text=paste0("sst.ECE$asst.",mmb,"=NULL")))
  
  }
}

#________________________________________
# Perform seasonal averages              \______________________________________________

# Define seasons
sst.ECE = sst.ECE[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
sst.ECE = sst.ECE[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
sst.ECE = sst.ECE[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
sst.ECE = sst.ECE[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

sst.ECE = sst.ECE[season == sel.season,]

# Add year for season
sst.ECE$s.year=year(sst.ECE$targetdate)
sst.ECE[targetmonth==12,]$s.year=sst.ECE[targetmonth==12,]$s.year+1


# Perform seasonal averages
for (mmb in members) {
  sst.ECE[ , eval(parse(text = paste0("s.asst.",mmb, ":=ave(asst.",mmb,")"))),by=.(s.year,lat,lon)]
}
sst.ECE = sst.ECE[,s.asst.em := ave(asst.em),by=.(s.year,lat,lon)]
# Eliminate monthly data
for (mmb in members) {
  eval(parse(text = paste0("sst.ECE$asst.",mmb, "=","NULL" )))
}
sst.ECE$asst.em=NULL

# Clean DT
sst.ECE$targetmonth=NULL
sst.ECE$season=NULL

# Remove repeated rows
sst.ECE = unique(sst.ECE, by=c("lat","lon","s.year"))


#________________________________________
# Perform ny-year rolling means           \______________________________________________

if(roll.years==TRUE){
  
  alig="left" #to be used for the rolling mean, using for time=t data between t and t+3
  
  sst.ECE=sst.ECE[, s.year.roll := frollmean(s.year,n=ny,align=alig),by=.(lat,lon)]
  sst.ECE=sst.ECE[, s.asst.em.roll := frollmean(s.asst.em,n=ny,align=alig),by=.(lat,lon)]
  for (mmb in members) {
    sst.ECE=sst.ECE[, eval(parse(text = paste0("s.asst.",mmb,".roll := frollmean(s.asst.",mmb,",n=ny,align=alig)"))),by=.(lat,lon)]
  }
  
  # Remove rows which could not be used to compute the rolling mean
  sst.ECE = sst.ECE[!is.na(s.year.roll),]
  # Eliminate previous variables and rearrange
  # sst.ECE$s.year=NULL
  # setnames(sst.ECE,"s.year.roll","s.year")
  sst.ECE$s.asst.em=NULL
  setnames(sst.ECE,"s.asst.em.roll","s.asst.em")
  for (mmb in members) {
    eval(parse(text = paste0("sst.ECE$s.asst.",mmb,"=NULL")))
    oldname=paste0("s.asst.",mmb,".roll")
    newname=paste0("s.asst.",mmb)
    
    eval(parse(text = paste0("setnames(sst.ECE,oldname,newname)")))
    
  }

}

#________________________________________
#  EOF calculation                       \______________________________________________

# Apply latitude weight to the anomalies
for (mmb in members) {
  eval(parse(text = paste0("sst.ECE$s.asst.",mmb, "=","sst.ECE$s.asst.",mmb,"*sqrt(cos(sst.ECE$lat*pi/180))" )))
}
sst.ECE$s.asst.em=sst.ECE$s.asst.em*sqrt(cos(sst.ECE$lat*pi/180))

# Create subset
sst.ECE.sub = sst.ECE

# Remove NA
sst.ECE.sub = sst.ECE.sub[!is.na(s.asst.em),]

# Rearrange a data table to use all 15 members to compute EOF
sst.ECE.allin1 = sst.ECE.sub[,c("lat","lon","targetdate")]
sst.ECE.2add = sst.ECE.allin1

for(reps in 1:14){
  sst.ECE.2add$targetdate=sst.ECE.2add$targetdate+(100*365) # Add 100 years to each round
  sst.ECE.allin1 = rbind(sst.ECE.allin1,sst.ECE.2add)
}

sst.ECE.allin1$asst=NA_real_

l=length(sst.ECE.sub$targetdate)

sst.ECE.allin1$asst[1:l]=sst.ECE.sub$s.asst.1
sst.ECE.allin1$asst[(1+l*1):(l*2)]=sst.ECE.sub$s.asst.2
sst.ECE.allin1$asst[(1+l*2):(l*3)]=sst.ECE.sub$s.asst.4
sst.ECE.allin1$asst[(1+l*3):(l*4)]=sst.ECE.sub$s.asst.5
sst.ECE.allin1$asst[(1+l*4):(l*5)]=sst.ECE.sub$s.asst.6
sst.ECE.allin1$asst[(1+l*5):(l*6)]=sst.ECE.sub$s.asst.7
sst.ECE.allin1$asst[(1+l*6):(l*7)]=sst.ECE.sub$s.asst.8
sst.ECE.allin1$asst[(1+l*7):(l*8)]=sst.ECE.sub$s.asst.9
sst.ECE.allin1$asst[(1+l*8):(l*9)]=sst.ECE.sub$s.asst.10
sst.ECE.allin1$asst[(1+l*9):(l*10)]=sst.ECE.sub$s.asst.11
sst.ECE.allin1$asst[(1+l*10):(l*11)]=sst.ECE.sub$s.asst.12
sst.ECE.allin1$asst[(1+l*11):(l*12)]=sst.ECE.sub$s.asst.13
sst.ECE.allin1$asst[(1+l*12):(l*13)]=sst.ECE.sub$s.asst.14
sst.ECE.allin1$asst[(1+l*13):(l*14)]=sst.ECE.sub$s.asst.15
sst.ECE.allin1$asst[(1+l*14):(l*15)]=sst.ECE.sub$s.asst.16


# Compute EOF for the ensemble mean 
# Computes EOFs 1 to 3 without extra latitude weight
eof = metR::EOF(asst ~ lat + lon | targetdate, n=1:3, data = sst.ECE.allin1, B = 100, probs = c(low = 0.1, hig = 0.9))

eof.sst.1 = cut(eof, 1)
eof.sst.2 = cut(eof, 2)
eof.sst.3 = cut(eof, 3)

# Normalize respect to standard deviation
eof.sst.1$right$pc1=eof.sst.1$right$asst/sd(eof.sst.1$right$asst,na.rm=TRUE)
eof.sst.2$right$pc2=eof.sst.2$right$asst/sd(eof.sst.2$right$asst,na.rm=TRUE)
eof.sst.3$right$pc3=eof.sst.3$right$asst/sd(eof.sst.3$right$asst,na.rm=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

sst.ECE.allin1=merge(sst.ECE.allin1,eof.sst.1$right[,.(targetdate,pc1)],by=c("targetdate"),all=TRUE)
sst.ECE.allin1=merge(sst.ECE.allin1,eof.sst.2$right[,.(targetdate,pc2)],by=c("targetdate"),all=TRUE)
sst.ECE.allin1=merge(sst.ECE.allin1,eof.sst.3$right[,.(targetdate,pc3)],by=c("targetdate"),all=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

sst.ECE.allin1=sst.ECE.allin1[,.(targetdate,asst,pc1,pc2,pc3,hcorrmap.eof1=cor(asst,pc1, use = "pairwise.complete.obs")),by=.(lat,lon)]
sst.ECE.allin1=sst.ECE.allin1[,.(targetdate,asst,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2=cor(asst,pc2, use = "pairwise.complete.obs")),by=.(lat,lon)]
sst.ECE.allin1=sst.ECE.allin1[,.(targetdate,asst,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3=cor(asst,pc3, use = "pairwise.complete.obs")),by=.(lat,lon)]

# Separate data tables for eofs and pcs
sst.eof=unique(sst.ECE.allin1[,.(lat,lon,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3)],by=c("lat","lon"))
sst.pcs.allin1=unique(sst.ECE.allin1[,.(targetdate,pc1,pc2,pc3)],by=c("targetdate"))

# I have 67 observations for each member, so dates and data table should be arrenged consequently
my=max(year(sst.ECE.sub$targetdate))
a=sum((year(sst.pcs.allin1$targetdate)==my)*seq(1,length(sst.pcs.allin1$targetdate),1))

sst.pcs=sst.pcs.allin1[1:a,]
setnames(sst.pcs,"pc1","pc1.1")
setnames(sst.pcs,"pc2","pc2.1")
setnames(sst.pcs,"pc3","pc3.1")

sst.pcs[,`:=` (pc1.2=sst.pcs.allin1$pc1[(1+a*1):(a*2)], pc2.2=sst.pcs.allin1$pc2[(1+a*1):(a*2)], pc3.2=sst.pcs.allin1$pc3[(1+a*1):(a*2)])]

for(i in c(4,5,6,7,8,9,10,11,12,13,14,15,16)){
  sst.pcs[, eval(parse(text=paste0("`:=` (pc1.",i,"=sst.pcs.allin1$pc1[(1+a*",(i-2),"):(a*",(i-1),")], pc2.",i,"=sst.pcs.allin1$pc2[(1+a*",(i-2),"):(a*",(i-1),")], pc3.",i,"=sst.pcs.allin1$pc3[(1+a*",(i-2),"):(a*",(i-1),")])"))) ]
}

sst.pcs$pc1.em = rowMeans(cbind(sst.pcs$pc1.1,sst.pcs$pc1.2,sst.pcs$pc1.4,sst.pcs$pc1.5,sst.pcs$pc1.6,sst.pcs$pc1.7,sst.pcs$pc1.8,sst.pcs$pc1.9,sst.pcs$pc1.10,sst.pcs$pc1.11,sst.pcs$pc1.12,sst.pcs$pc1.13,sst.pcs$pc1.14,sst.pcs$pc1.15,sst.pcs$pc1.16))
sst.pcs$pc2.em = rowMeans(cbind(sst.pcs$pc2.1,sst.pcs$pc2.2,sst.pcs$pc2.4,sst.pcs$pc2.5,sst.pcs$pc2.6,sst.pcs$pc2.7,sst.pcs$pc2.8,sst.pcs$pc2.9,sst.pcs$pc2.10,sst.pcs$pc2.11,sst.pcs$pc2.12,sst.pcs$pc2.13,sst.pcs$pc2.14,sst.pcs$pc2.15,sst.pcs$pc2.16))
sst.pcs$pc3.em = rowMeans(cbind(sst.pcs$pc3.1,sst.pcs$pc3.2,sst.pcs$pc3.4,sst.pcs$pc3.5,sst.pcs$pc3.6,sst.pcs$pc3.7,sst.pcs$pc3.8,sst.pcs$pc3.9,sst.pcs$pc3.10,sst.pcs$pc3.11,sst.pcs$pc3.12,sst.pcs$pc3.13,sst.pcs$pc3.14,sst.pcs$pc3.15,sst.pcs$pc3.16))


  
#________________________________________
# Plotting and saving                    \______________________________________________

# Change longitude to -180:180 to plot correctly
sst.eof[lon>180]$lon=sst.eof[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-0.6
bmax=0.6
bstep=0.1
bbreaks=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6)

g1 = ggplot() +
  geom_contour_fill(data=sst.eof,aes(lon, lat, z = hcorrmap.eof1),na.fill=TRUE)+
  scale_fill_distiller(name="EOF1",palette="RdBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF1. Explained variance: ",as.character(round(eof.sst.1$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  geom_contour_fill(data=sst.eof,aes(lon, lat, z = hcorrmap.eof2),na.fill=TRUE)+
  scale_fill_distiller(name="EOF2",palette="RdBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF2. Explained variance: ",as.character(round(eof.sst.2$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g3 = ggplot() +
  geom_contour_fill(data=sst.eof,aes(lon, lat, z = hcorrmap.eof3),na.fill=TRUE)+
  scale_fill_distiller(name="EOF3",palette="RdBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF3. Explained variance: ",as.character(round(eof.sst.3$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))



g4 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(targetdate, pc1.1),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.2),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.4),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.5),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.6),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.7),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.8),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.9),color="#969696",alpha=0.6)+  
  geom_line(data=sst.pcs,aes(targetdate, pc1.10),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.11),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.12),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.13),color="#969696",alpha=0.6)+  
  geom_line(data=sst.pcs,aes(targetdate, pc1.14),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.15),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.16),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc1.em),color="#dd1c77",size=1.25)+
  scale_x_date(breaks=sst.pcs$targetdate[seq(1,length(sst.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC1 (normalized units)")+
  theme(axis.title = element_text(size=11))


g5 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(targetdate, pc2.1),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.2),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.4),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.5),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.6),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.7),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.8),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.9),color="#969696",alpha=0.6)+  
  geom_line(data=sst.pcs,aes(targetdate, pc2.10),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.11),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.12),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.13),color="#969696",alpha=0.6)+  
  geom_line(data=sst.pcs,aes(targetdate, pc2.14),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.15),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.16),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc2.em),color="#dd1c77",size=1.25)+
  scale_x_date(breaks=sst.pcs$targetdate[seq(1,length(sst.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC2 (normalized units)")+
  theme(axis.title = element_text(size=11))

g6 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(targetdate, pc3.1),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.2),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.4),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.5),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.6),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.7),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.8),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.9),color="#969696",alpha=0.6)+  
  geom_line(data=sst.pcs,aes(targetdate, pc3.10),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.11),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.12),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.13),color="#969696",alpha=0.6)+  
  geom_line(data=sst.pcs,aes(targetdate, pc3.14),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.15),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.16),color="#969696",alpha=0.6)+
  geom_line(data=sst.pcs,aes(targetdate, pc3.em),color="#dd1c77",size=1.25)+
  scale_x_date(breaks=sst.pcs$targetdate[seq(1,length(sst.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC3 (normalized units)")+
  theme(axis.title = element_text(size=11))

# Save Figure and PCs for the lead

if(remove.trend==TRUE){
  # save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/",sel.season,"_historical_ECEarth3_PCs_EOFs_SST_weighted_notrend_allmembers.RData"))
  if(roll.years==TRUE){
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of SST anomalies (no trend, rolling ",ny," years, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_SST_weighted_notrend_allmembers_rollingyears_",ny,".png"),plot=fig,width = 10, height = 8)
    
  }else{
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of SST anomalies (no trend, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_SST_weighted_notrend_allmembers.png"),plot=fig,width = 10, height = 8)
    
  }
  
  }else{
  # save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/",sel.season,"_historical_ECEarth3_PCs_EOFs_SST_weighted_allmembers.RData"))
  if(roll.years==TRUE){
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of SST anomalies (rolling ",ny," years, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_SST_weighted_allmembers_rollingyears_",ny,".png"),plot=fig,width = 10, height = 8)
    
  }else{
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of SST anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_SST_weighted_allmembers.png"),plot=fig,width = 10, height = 8)
    
  }
} 

