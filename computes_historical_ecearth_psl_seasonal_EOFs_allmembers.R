# This script opens the initialized decadal predictions of PSL (tos) of 
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

#---------------------------------------------------------------------------------------
#  Man Program
#---------------------------------------------------------------------------------------

#_________________________________________
# Loads data                              \_____________________________________________

load(file="/home/maralv/data/apsl.ECEarth.hist.19492014.RData")

#################### Settings ########################
# Remove trend? TRUE or FALSE
remove.trend=FALSE
# Select season
sel.season="DJF"
# Rolling years
roll.years=FALSE
# number of years to roll
ny=4
# Members being used
members=c(1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)
######################################################

#________________________________________
# Linear trend removal while monthly     \______________________________________________

# For now, only for ensemble mean respect to itself

if(remove.trend==TRUE){
  # Remove linear trend of PSL anomalies
  psl.ECE=psl.ECE[order(targetdate),]
  psl.ECE=psl.ECE[,dt.apsl.em := detrend(apsl.em),by=.(lat,lon)]
  psl.ECE$apsl.em=NULL
  setnames(psl.ECE,"dt.apsl.em","apsl.em")
  
  for (mmb in members) {
  psl.ECE[ , eval(parse(text = paste0("apsl.",mmb, ":= detrend(apsl.",mmb,")"))),by=.(lat,lon)]
  # eval(parse(text=paste0("psl.ECE$apsl.",mmb,"=NULL")))
  
  }
}

#________________________________________
# Perform seasonal averages              \______________________________________________

# Define seasons
psl.ECE = psl.ECE[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
psl.ECE = psl.ECE[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
psl.ECE = psl.ECE[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
psl.ECE = psl.ECE[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

psl.ECE = psl.ECE[season == sel.season,]

# Add year for season
psl.ECE$s.year=year(psl.ECE$targetdate)
psl.ECE[targetmonth==12,]$s.year=psl.ECE[targetmonth==12,]$s.year+1


# Perform seasonal averages
for (mmb in members) {
  psl.ECE[ , eval(parse(text = paste0("s.apsl.",mmb, ":=ave(apsl.",mmb,")"))),by=.(s.year,lat,lon)]
}
psl.ECE = psl.ECE[,s.apsl.em := ave(apsl.em),by=.(s.year,lat,lon)]
# Eliminate monthly data
for (mmb in members) {
  eval(parse(text = paste0("psl.ECE$apsl.",mmb, "=","NULL" )))
}
psl.ECE$apsl.em=NULL

# Clean DT
psl.ECE$targetmonth=NULL
psl.ECE$season=NULL

# Remove repeated rows
psl.ECE = unique(psl.ECE, by=c("lat","lon","s.year"))


#________________________________________
# Perform ny-year rolling means           \______________________________________________

if(roll.years==TRUE){
  
  alig="left" #to be used for the rolling mean, using for time=t data between t and t+3
  
  psl.ECE=psl.ECE[, s.year.roll := frollmean(s.year,n=ny,align=alig),by=.(lat,lon)]
  psl.ECE=psl.ECE[, s.apsl.em.roll := frollmean(s.apsl.em,n=ny,align=alig),by=.(lat,lon)]
  for (mmb in members) {
    psl.ECE=psl.ECE[, eval(parse(text = paste0("s.apsl.",mmb,".roll := frollmean(s.apsl.",mmb,",n=ny,align=alig)"))),by=.(lat,lon)]
  }
  
  # Remove rows which could not be used to compute the rolling mean
  psl.ECE = psl.ECE[!is.na(s.year.roll),]
  # Eliminate previous variables and rearrange
  # psl.ECE$s.year=NULL
  # setnames(psl.ECE,"s.year.roll","s.year")
  psl.ECE$s.apsl.em=NULL
  setnames(psl.ECE,"s.apsl.em.roll","s.apsl.em")
  for (mmb in members) {
    eval(parse(text = paste0("psl.ECE$s.apsl.",mmb,"=NULL")))
    oldname=paste0("s.apsl.",mmb,".roll")
    newname=paste0("s.apsl.",mmb)
    
    eval(parse(text = paste0("setnames(psl.ECE,oldname,newname)")))
    
  }

}

#________________________________________
#  EOF calculation                       \______________________________________________

# Apply latitude weight to the anomalies
for (mmb in members) {
  eval(parse(text = paste0("psl.ECE$s.apsl.",mmb, "=","psl.ECE$s.apsl.",mmb,"*sqrt(cos(psl.ECE$lat*pi/180))" )))
}
psl.ECE$s.apsl.em=psl.ECE$s.apsl.em*sqrt(cos(psl.ECE$lat*pi/180))

# Create subset
psl.ECE.sub = psl.ECE

# Remove NA
psl.ECE.sub = psl.ECE.sub[!is.na(s.apsl.em),]

# Rearrange a data table to use all 15 members to compute EOF
psl.ECE.allin1 = psl.ECE.sub[,c("lat","lon","targetdate")]
psl.ECE.2add = psl.ECE.allin1

for(reps in 1:14){
  psl.ECE.2add$targetdate=psl.ECE.2add$targetdate+(100*365) # Add 100 years to each round
  psl.ECE.allin1 = rbind(psl.ECE.allin1,psl.ECE.2add)
}

psl.ECE.allin1$apsl=NA_real_

l=length(psl.ECE.sub$targetdate)

psl.ECE.allin1$apsl[1:l]=psl.ECE.sub$s.apsl.1
psl.ECE.allin1$apsl[(1+l*1):(l*2)]=psl.ECE.sub$s.apsl.2
psl.ECE.allin1$apsl[(1+l*2):(l*3)]=psl.ECE.sub$s.apsl.4
psl.ECE.allin1$apsl[(1+l*3):(l*4)]=psl.ECE.sub$s.apsl.5
psl.ECE.allin1$apsl[(1+l*4):(l*5)]=psl.ECE.sub$s.apsl.6
psl.ECE.allin1$apsl[(1+l*5):(l*6)]=psl.ECE.sub$s.apsl.7
psl.ECE.allin1$apsl[(1+l*6):(l*7)]=psl.ECE.sub$s.apsl.8
psl.ECE.allin1$apsl[(1+l*7):(l*8)]=psl.ECE.sub$s.apsl.9
psl.ECE.allin1$apsl[(1+l*8):(l*9)]=psl.ECE.sub$s.apsl.10
psl.ECE.allin1$apsl[(1+l*9):(l*10)]=psl.ECE.sub$s.apsl.11
psl.ECE.allin1$apsl[(1+l*10):(l*11)]=psl.ECE.sub$s.apsl.12
psl.ECE.allin1$apsl[(1+l*11):(l*12)]=psl.ECE.sub$s.apsl.13
psl.ECE.allin1$apsl[(1+l*12):(l*13)]=psl.ECE.sub$s.apsl.14
psl.ECE.allin1$apsl[(1+l*13):(l*14)]=psl.ECE.sub$s.apsl.15
psl.ECE.allin1$apsl[(1+l*14):(l*15)]=psl.ECE.sub$s.apsl.16


# Compute EOF for the ensemble mean 
# Computes EOFs 1 to 3 without extra latitude weight
eof = metR::EOF(apsl ~ lat + lon | targetdate, n=1:3, data = psl.ECE.allin1, B = 100, probs = c(low = 0.1, hig = 0.9))

eof.psl.1 = cut(eof, 1)
eof.psl.2 = cut(eof, 2)
eof.psl.3 = cut(eof, 3)

# Normalize respect to standard deviation
eof.psl.1$right$pc1=eof.psl.1$right$apsl/sd(eof.psl.1$right$apsl,na.rm=TRUE)
eof.psl.2$right$pc2=eof.psl.2$right$apsl/sd(eof.psl.2$right$apsl,na.rm=TRUE)
eof.psl.3$right$pc3=eof.psl.3$right$apsl/sd(eof.psl.3$right$apsl,na.rm=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the psl anomaly at that grid point.

psl.ECE.allin1=merge(psl.ECE.allin1,eof.psl.1$right[,.(targetdate,pc1)],by=c("targetdate"),all=TRUE)
psl.ECE.allin1=merge(psl.ECE.allin1,eof.psl.2$right[,.(targetdate,pc2)],by=c("targetdate"),all=TRUE)
psl.ECE.allin1=merge(psl.ECE.allin1,eof.psl.3$right[,.(targetdate,pc3)],by=c("targetdate"),all=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the psl anomaly at that grid point.

psl.ECE.allin1=psl.ECE.allin1[,.(targetdate,apsl,pc1,pc2,pc3,hcorrmap.eof1=cor(apsl,pc1, use = "pairwise.complete.obs")),by=.(lat,lon)]
psl.ECE.allin1=psl.ECE.allin1[,.(targetdate,apsl,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2=cor(apsl,pc2, use = "pairwise.complete.obs")),by=.(lat,lon)]
psl.ECE.allin1=psl.ECE.allin1[,.(targetdate,apsl,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3=cor(apsl,pc3, use = "pairwise.complete.obs")),by=.(lat,lon)]

# Separate data tables for eofs and pcs
psl.eof=unique(psl.ECE.allin1[,.(lat,lon,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3)],by=c("lat","lon"))
psl.pcs.allin1=unique(psl.ECE.allin1[,.(targetdate,pc1,pc2,pc3)],by=c("targetdate"))

# I have 67 observations for each member, so dates and data table should be arrenged consequently
my=max(year(psl.ECE.sub$targetdate))
a=sum((year(psl.pcs.allin1$targetdate)==my)*seq(1,length(psl.pcs.allin1$targetdate),1))

psl.pcs=psl.pcs.allin1[1:a,]
setnames(psl.pcs,"pc1","pc1.1")
setnames(psl.pcs,"pc2","pc2.1")
setnames(psl.pcs,"pc3","pc3.1")

psl.pcs[,`:=` (pc1.2=psl.pcs.allin1$pc1[(1+a*1):(a*2)], pc2.2=psl.pcs.allin1$pc2[(1+a*1):(a*2)], pc3.2=psl.pcs.allin1$pc3[(1+a*1):(a*2)])]

for(i in c(4,5,6,7,8,9,10,11,12,13,14,15,16)){
  psl.pcs[, eval(parse(text=paste0("`:=` (pc1.",i,"=psl.pcs.allin1$pc1[(1+a*",(i-2),"):(a*",(i-1),")], pc2.",i,"=psl.pcs.allin1$pc2[(1+a*",(i-2),"):(a*",(i-1),")], pc3.",i,"=psl.pcs.allin1$pc3[(1+a*",(i-2),"):(a*",(i-1),")])"))) ]
}

psl.pcs$pc1.em = rowMeans(cbind(psl.pcs$pc1.1,psl.pcs$pc1.2,psl.pcs$pc1.4,psl.pcs$pc1.5,psl.pcs$pc1.6,psl.pcs$pc1.7,psl.pcs$pc1.8,psl.pcs$pc1.9,psl.pcs$pc1.10,psl.pcs$pc1.11,psl.pcs$pc1.12,psl.pcs$pc1.13,psl.pcs$pc1.14,psl.pcs$pc1.15,psl.pcs$pc1.16))
psl.pcs$pc2.em = rowMeans(cbind(psl.pcs$pc2.1,psl.pcs$pc2.2,psl.pcs$pc2.4,psl.pcs$pc2.5,psl.pcs$pc2.6,psl.pcs$pc2.7,psl.pcs$pc2.8,psl.pcs$pc2.9,psl.pcs$pc2.10,psl.pcs$pc2.11,psl.pcs$pc2.12,psl.pcs$pc2.13,psl.pcs$pc2.14,psl.pcs$pc2.15,psl.pcs$pc2.16))
psl.pcs$pc3.em = rowMeans(cbind(psl.pcs$pc3.1,psl.pcs$pc3.2,psl.pcs$pc3.4,psl.pcs$pc3.5,psl.pcs$pc3.6,psl.pcs$pc3.7,psl.pcs$pc3.8,psl.pcs$pc3.9,psl.pcs$pc3.10,psl.pcs$pc3.11,psl.pcs$pc3.12,psl.pcs$pc3.13,psl.pcs$pc3.14,psl.pcs$pc3.15,psl.pcs$pc3.16))


  
#________________________________________
# Plotting and saving                    \______________________________________________

# Change longitude to -180:180 to plot correctly
psl.eof[lon>180]$lon=psl.eof[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-0.6
bmax=0.6
bstep=0.1
bbreaks.contours=c(-99,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,99)
bbreaks.cbar=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6)
labels.cbar=as.character(bbreaks.cbar)
labels.cbar[1]=""
labels.cbar[length(labels.cbar)]=""
  
g1 = ggplot() +
  geom_contour_fill(data=psl.eof,aes(lon, lat, z = hcorrmap.eof1),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF1",palette="RdYlBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF1. Explained variance: ",as.character(round(eof.psl.1$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  geom_contour_fill(data=psl.eof,aes(lon, lat, z = hcorrmap.eof2),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF2",palette="RdYlBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF2. Explained variance: ",as.character(round(eof.psl.2$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g3 = ggplot() +
  geom_contour_fill(data=psl.eof,aes(lon, lat, z = hcorrmap.eof3),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="EOF3",palette="RdYlBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF3. Explained variance: ",as.character(round(eof.psl.3$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))



g4 = ggplot() +
  theme_bw()+
  geom_line(data=psl.pcs,aes(targetdate, pc1.1),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.2),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.4),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.5),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.6),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.7),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.8),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.9),color="#969696",alpha=0.6)+  
  geom_line(data=psl.pcs,aes(targetdate, pc1.10),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.11),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.12),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.13),color="#969696",alpha=0.6)+  
  geom_line(data=psl.pcs,aes(targetdate, pc1.14),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.15),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.16),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc1.em),color="#dd1c77",size=1.25)+
  scale_x_date(breaks=psl.pcs$targetdate[seq(1,length(psl.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC1 (normalized units)")+
  theme(axis.title = element_text(size=11))


g5 = ggplot() +
  theme_bw()+
  geom_line(data=psl.pcs,aes(targetdate, pc2.1),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.2),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.4),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.5),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.6),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.7),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.8),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.9),color="#969696",alpha=0.6)+  
  geom_line(data=psl.pcs,aes(targetdate, pc2.10),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.11),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.12),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.13),color="#969696",alpha=0.6)+  
  geom_line(data=psl.pcs,aes(targetdate, pc2.14),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.15),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.16),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc2.em),color="#dd1c77",size=1.25)+
  scale_x_date(breaks=psl.pcs$targetdate[seq(1,length(psl.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC2 (normalized units)")+
  theme(axis.title = element_text(size=11))

g6 = ggplot() +
  theme_bw()+
  geom_line(data=psl.pcs,aes(targetdate, pc3.1),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.2),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.4),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.5),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.6),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.7),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.8),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.9),color="#969696",alpha=0.6)+  
  geom_line(data=psl.pcs,aes(targetdate, pc3.10),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.11),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.12),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.13),color="#969696",alpha=0.6)+  
  geom_line(data=psl.pcs,aes(targetdate, pc3.14),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.15),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.16),color="#969696",alpha=0.6)+
  geom_line(data=psl.pcs,aes(targetdate, pc3.em),color="#dd1c77",size=1.25)+
  scale_x_date(breaks=psl.pcs$targetdate[seq(1,length(psl.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC3 (normalized units)")+
  theme(axis.title = element_text(size=11))

# Save Figure and PCs for the lead

if(remove.trend==TRUE){
  # save(psl.pcs,psl.eof,file=paste0("/home/maralv/data/",sel.season,"_historical_ECEarth3_PCs_EOFs_PSL_weighted_notrend_allmembers.RData"))
  if(roll.years==TRUE){
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of PSL anomalies (no trend, rolling ",ny," years, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_PSL_weighted_notrend_allmembers_rollingyears_",ny,".png"),plot=fig,width = 10, height = 8)
    
  }else{
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of PSL anomalies (no trend, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_PSL_weighted_notrend_allmembers.png"),plot=fig,width = 10, height = 8)
    
  }
  
  }else{
  # save(psl.pcs,psl.eof,file=paste0("/home/maralv/data/",sel.season,"_historical_ECEarth3_PCs_EOFs_PSL_weighted_allmembers.RData"))
  if(roll.years==TRUE){
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of PSL anomalies (rolling ",ny," years, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_PSL_weighted_allmembers_rollingyears_",ny,".png"),plot=fig,width = 10, height = 8)
    
  }else{
    fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , historical: EOFs of PSL anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_historical_EOFs_PSL_weighted_allmembers.png"),plot=fig,width = 10, height = 8)
    
  }
} 

