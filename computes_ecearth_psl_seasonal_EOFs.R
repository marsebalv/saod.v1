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

load(file="/home/maralv/data/apsl.ECEarth.DCP.19612017.RData")

#________________________________________
# Perform seasonal averages              \______________________________________________

# Select season for computation:
sel.season = "DJF"

# Define seasons
psl.ECE = psl.ECE[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
psl.ECE = psl.ECE[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
psl.ECE = psl.ECE[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
psl.ECE = psl.ECE[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

#!!!! BEWARE NOT TO USE LEAD.YEAR=0 WHEN USING SON!!!!!

psl.ECE = psl.ECE[season == sel.season,]

# Perform seasonal averages

for (mmb in 1:15) {
  psl.ECE[ , eval(parse(text = paste0("s.apsl.",mmb, ":=ave(apsl.",mmb,")"))),by=.(lead.year,startdate,lat,lon)]
}
psl.ECE = psl.ECE[,s.apsl.em := ave(apsl.em),by=.(lead.year,startdate,lat,lon)]
# Eliminate monthly data
for (mmb in 1:15) {
  eval(parse(text = paste0("psl.ECE$apsl.",mmb, "=","NULL" )))
}
psl.ECE$apsl.em=NULL

# Clean DT
psl.ECE$targetmonth=NULL
psl.ECE$lead=NULL
psl.ECE$season=NULL

# Remove repeated rows
psl.ECE = unique(psl.ECE, by=c("lat","lon","lead.year","startdate"))

psl.ECE = psl.ECE[lead.year<11,]

#________________________________________
# Lead selection and EOF calculation     \______________________________________________

sel.lead=7

# Apply latitude weight to the anomalies
for (mmb in 1:15) {
  eval(parse(text = paste0("psl.ECE$s.apsl.",mmb, "=","psl.ECE$s.apsl.",mmb,"*sqrt(cos(psl.ECE$lat*pi/180))" )))
}
psl.ECE$s.apsl.em=psl.ECE$s.apsl.em*sqrt(cos(psl.ECE$lat*pi/180))

# Select lead to compute EOF
psl.ECE.sub = psl.ECE[lead.year==sel.lead,]

# Compute EOF for the ensemble mean psl
# Computes EOFs 1 to 3 without latitude weight
eof.slp = metR::EOF(s.apsl.em ~ lat + lon | targetdate, n=1:3, data = psl.ECE.sub, B = 100, probs = c(low = 0.1, hig = 0.9))

eof.slp.1 = cut(eof.slp, 1)
eof.slp.2 = cut(eof.slp, 2)
eof.slp.3 = cut(eof.slp, 3)

# Normalize respect to standard deviation
eof.slp.1$right$pc1=eof.slp.1$right$s.apsl.em/sd(eof.slp.1$right$s.apsl.em,na.rm=TRUE)
eof.slp.2$right$pc2=eof.slp.2$right$s.apsl.em/sd(eof.slp.2$right$s.apsl.em,na.rm=TRUE)
eof.slp.3$right$pc3=eof.slp.3$right$s.apsl.em/sd(eof.slp.3$right$s.apsl.em,na.rm=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

psl.ECE.sub=merge(psl.ECE.sub,eof.slp.1$right[,.(targetdate,pc1)],by=c("targetdate"),all=TRUE)
psl.ECE.sub=merge(psl.ECE.sub,eof.slp.2$right[,.(targetdate,pc2)],by=c("targetdate"),all=TRUE)
psl.ECE.sub=merge(psl.ECE.sub,eof.slp.3$right[,.(targetdate,pc3)],by=c("targetdate"),all=TRUE)

# For the spatial patterns: presented as homogeneous correlation maps; i.e., the contours are scaled such
# that the value at each grid point is the correlation coefficient between the time series of expansion
# coefficients of Fig. 5 and the SST anomaly at that grid point.

psl.ECE.sub=psl.ECE.sub[,.(targetdate,s.apsl.em,pc1,pc2,pc3,hcorrmap.eof1=cor(s.apsl.em,pc1, use = "pairwise.complete.obs")),by=.(lat,lon)]
psl.ECE.sub=psl.ECE.sub[,.(targetdate,s.apsl.em,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2=cor(s.apsl.em,pc2, use = "pairwise.complete.obs")),by=.(lat,lon)]
psl.ECE.sub=psl.ECE.sub[,.(targetdate,s.apsl.em,pc1,pc2,pc3,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3=cor(s.apsl.em,pc3, use = "pairwise.complete.obs")),by=.(lat,lon)]

# Separate data tables for eofs and pcs
slp.eof=unique(psl.ECE.sub[,.(lat,lon,hcorrmap.eof1,hcorrmap.eof2,hcorrmap.eof3)],by=c("lat","lon"))
slp.pcs=unique(psl.ECE.sub[,.(targetdate,pc1,pc2,pc3)],by=c("targetdate"))


#________________________________________
# Plotting and saving                    \______________________________________________

# Change longitude to -180:180 to plot correctly
slp.eof[lon>180]$lon=slp.eof[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-1
bmax=0.5
bstep=0.05
#bbreaks=seq(-1,0.1,0.05)
bbreaks=c(-1,-0.95,-0.9,-0.85,-0.8,-0.7,-0.5,-0.2,0,0.2,0.3,0.4,0.5)

g1 = ggplot() +
  geom_contour_fill(data=slp.eof,aes(lon, lat, z = hcorrmap.eof1),na.fill=TRUE)+
  scale_fill_distiller(name="EOF1",palette="RdYlBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF1. Explained variance: ",as.character(round(eof.slp.1$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  geom_contour_fill(data=slp.eof,aes(lon, lat, z = hcorrmap.eof2),na.fill=TRUE)+
  scale_fill_distiller(name="EOF2",palette="RdYlBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF2. Explained variance: ",as.character(round(eof.slp.2$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g3 = ggplot() +
  geom_contour_fill(data=slp.eof,aes(lon, lat, z = hcorrmap.eof3),na.fill=TRUE)+
  scale_fill_distiller(name="EOF3",palette="RdYlBu",direction=-1,
                       breaks=bbreaks,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       oob  = scales::squish)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("EOF3. Explained variance: ",as.character(round(eof.slp.3$sdev$r2*100,1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))



g4 = ggplot() +
  theme_bw()+
  geom_line(data=slp.pcs,aes(targetdate, pc1))+
  scale_x_date(breaks=slp.pcs$targetdate[seq(1,length(slp.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC1 (normalized units)")+
  theme(axis.title = element_text(size=11))

g5 = ggplot() +
  theme_bw()+
  geom_line(data=slp.pcs,aes(targetdate, pc2))+
  scale_x_date(breaks=slp.pcs$targetdate[seq(1,length(slp.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC2 (normalized units)")+
  theme(axis.title = element_text(size=11))

g6 = ggplot() +
  theme_bw()+
  geom_line(data=slp.pcs,aes(targetdate, pc3))+
  scale_x_date(breaks=slp.pcs$targetdate[seq(1,length(slp.pcs$targetdate),5)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3.75,3.75),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("PC3 (normalized units)")+
  theme(axis.title = element_text(size=11))



fig <- grid.arrange(g1,g4,g2,g5,g3,g6, ncol = 2,top = textGrob(paste0(sel.season," , lead = ",sel.lead,": EOFs of SLP anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"lead",sel.lead,"_EOFs_SLPA_weighted.png"),plot=fig,width = 10, height = 8)