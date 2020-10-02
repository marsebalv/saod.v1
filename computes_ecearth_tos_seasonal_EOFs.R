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

load(file="/home/maralv/data/asst.ECEarth.DCP.19612017.RData")

#________________________________________
# Perform seasonal averages              \______________________________________________

# Select season for computation:
sel.season = "DJF"

# Define seasons
sst.ECE = sst.ECE[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
sst.ECE = sst.ECE[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
sst.ECE = sst.ECE[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
sst.ECE = sst.ECE[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

#!!!! BEWARE NOT TO USE LEAD.YEAR=0 WHEN USING SON!!!!!

sst.ECE = sst.ECE[season == sel.season,]

# Perform seasonal averages

for (mmb in 1:15) {
  sst.ECE[ , eval(parse(text = paste0("s.asst.",mmb, ":=ave(asst.",mmb,")"))),by=.(lead.year,startdate,lat,lon)]
}
sst.ECE = sst.ECE[,s.asst.em := ave(asst.em),by=.(lead.year,startdate,lat,lon)]
# Eliminate monthly data
for (mmb in 1:15) {
  eval(parse(text = paste0("sst.ECE$asst.",mmb, "=","NULL" )))
}
sst.ECE$asst.em=NULL

# Clean DT
sst.ECE$targetmonth=NULL
sst.ECE$lead=NULL
sst.ECE$season=NULL

# Remove repeated rows
sst.ECE = unique(sst.ECE, by=c("lat","lon","lead.year","startdate"))

sst.ECE = sst.ECE[lead.year<11,]

#________________________________________
# Lead selection and EOF calculation     \______________________________________________

sel.lead=1

# Apply latitude weight to the anomalies
for (mmb in 1:15) {
  eval(parse(text = paste0("sst.ECE$s.asst.",mmb, "=","sst.ECE$s.asst.",mmb,"*sqrt(cos(sst.ECE$lat*pi/180))" )))
}
sst.ECE$s.asst.em=sst.ECE$s.asst.em*sqrt(cos(sst.ECE$lat*pi/180))

# Select lead to compute EOF
sst.ECE.sub = sst.ECE[lead.year==sel.lead,]

# Remove NA
sst.ECE.sub = sst.ECE.sub[!is.na(s.asst.em),]

# Compute EOF for the ensemble mean psl
# Computes EOFs 1 to 3 without latitude weight
eof.slp = metR::EOF(s.asst.em ~ lat + lon | targetdate, n=1:3, data = sst.ECE.sub, B = 100, probs = c(low = 0.1, hig = 0.9))

eof.sst.1 = cut(eof.slp, 1)
eof.sst.2 = cut(eof.slp, 2)
eof.sst.3 = cut(eof.slp, 3)

# Normalize respect to standard deviation
eof.sst.1$right$pc1=eof.sst.1$right$s.asst.em/sd(eof.sst.1$right$s.asst.em,na.rm=TRUE)
eof.sst.2$right$pc2=eof.sst.2$right$s.asst.em/sd(eof.sst.2$right$s.asst.em,na.rm=TRUE)
eof.sst.3$right$pc3=eof.sst.3$right$s.asst.em/sd(eof.sst.3$right$s.asst.em,na.rm=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

sst.ECE.sub=merge(sst.ECE.sub,eof.sst.1$right[,.(targetdate,pc1)],by=c("targetdate"),all=TRUE)
sst.ECE.sub=merge(sst.ECE.sub,eof.sst.2$right[,.(targetdate,pc2)],by=c("targetdate"),all=TRUE)
sst.ECE.sub=merge(sst.ECE.sub,eof.sst.3$right[,.(targetdate,pc3)],by=c("targetdate"),all=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

sst.ECE.sub=sst.ECE.sub[,.(targetdate,s.asst.em,pc1,pc2,pc3,hcorrmap.eof1=cor(s.asst.em,pc1, use = "pairwise.complete.obs")),by=.(lat,lon)]
sst.ECE.sub=sst.ECE.sub[,.(targetdate,s.asst.em,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2=cor(s.asst.em,pc2, use = "pairwise.complete.obs")),by=.(lat,lon)]
sst.ECE.sub=sst.ECE.sub[,.(targetdate,s.asst.em,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3=cor(s.asst.em,pc3, use = "pairwise.complete.obs")),by=.(lat,lon)]

# Separate data tables for eofs and pcs
sst.eof=unique(sst.ECE.sub[,.(lat,lon,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3)],by=c("lat","lon"))
sst.pcs=unique(sst.ECE.sub[,.(targetdate,pc1,pc2,pc3)],by=c("targetdate"))


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
  geom_line(data=sst.pcs,aes(targetdate, pc1))+
  scale_x_date(breaks=sst.pcs$targetdate[seq(1,length(sst.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC1 (normalized units)")+
  theme(axis.title = element_text(size=11))

g5 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(targetdate, pc2))+
  scale_x_date(breaks=sst.pcs$targetdate[seq(1,length(sst.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC2 (normalized units)")+
  theme(axis.title = element_text(size=11))

g6 = ggplot() +
  theme_bw()+
  geom_line(data=sst.pcs,aes(targetdate, pc3))+
  scale_x_date(breaks=sst.pcs$targetdate[seq(1,length(sst.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC3 (normalized units)")+
  theme(axis.title = element_text(size=11))


fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , lead = ",sel.lead,": EOFs of SST anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"lead",sel.lead,"_EOFs_SST_weighted.png"),plot=fig,width = 10, height = 8)

