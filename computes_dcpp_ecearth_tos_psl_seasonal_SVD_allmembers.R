# This script opens the initialized decadal predictions of SST (tos) & SLP (psl) of 
# EC-Earth CMIP6 and computes the SAOD pattern (SVD1)
#
# M. Alvarez (2021)
#--------------------------------------------------------------------------------

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

rm(list=ls())
graphics.off()

# Start for on leads
for(l in 1:10){

setwd("/home/maralv/")

#---------------------------------------------------------------------------------------
#  functions 
#---------------------------------------------------------------------------------------

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



#################### Settings ########################
# # Remove trend? TRUE or FALSE
remove.trend=FALSE
# Select season
sel.season="DJF"
# Select lead
sel.lead=l
######################################################

#_________________________________________
# Loads data                              \_____________________________________________

load(file="/home/maralv/data/apsl.ECEarth.DCP.19612017.RData")
load(file="/home/maralv/data/asst.ECEarth.DCP.19612017.RData")

# Select seasons first to reduce the dataset
# Define seasons
# SST
sst.ECE = sst.ECE[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
sst.ECE = sst.ECE[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
sst.ECE = sst.ECE[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
sst.ECE = sst.ECE[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

sst.ECE = sst.ECE[season == sel.season,]
sst.ECE$season = NULL # Not needed anymore

# SLP
psl.ECE = psl.ECE[targetmonth==12 | targetmonth==1 | targetmonth==2, season := "DJF" ]
psl.ECE = psl.ECE[targetmonth==3 | targetmonth==4 | targetmonth==5, season := "MAM" ]
psl.ECE = psl.ECE[targetmonth==6 | targetmonth==7 | targetmonth==8, season := "JJA" ]
psl.ECE = psl.ECE[targetmonth==9 | targetmonth==10 | targetmonth==11, season := "SON" ]

psl.ECE = psl.ECE[season == sel.season,]
psl.ECE$season = NULL # Not needed anymore

# Creo que ya "lead" no lo necesitaria porque usaria lead.year, y puedo tambien seleccionar ahora el lead con el que voy a trabajar.
# Reduce dataset
# SST
sst.ECE$lead = NULL # No more month detail needed
sst.ECE = sst.ECE[lead.year == sel.lead,]
sst.ECE$lead.year = NULL
sst.ECE$targetmonth = NULL

# SLP
psl.ECE$lead = NULL # No more month detail needed
psl.ECE = psl.ECE[lead.year == sel.lead,]
psl.ECE$lead.year = NULL
psl.ECE$targetmonth = NULL


# Rearrange data using number of member as variable
# STT
# rename according to member
old = c("asst.1","asst.2","asst.3","asst.4","asst.5","asst.6","asst.7","asst.8","asst.9","asst.10","asst.11","asst.12","asst.13","asst.14","asst.15","asst.em")
new = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","99")
setnames(sst.ECE,old,new)
# Melt data table for easy operation
sst.ECE=melt(sst.ECE,id=c("lat","lon","startdate","targetdate"),variable.name="member",value.name="asst")

sst.ECE$member=as.numeric(as.character(sst.ECE$member))

# SLP
# rename according to member
old = c("apsl.1","apsl.2","apsl.3","apsl.4","apsl.5","apsl.6","apsl.7","apsl.8","apsl.9","apsl.10","apsl.11","apsl.12","apsl.13","apsl.14","apsl.15","apsl.em")
setnames(psl.ECE,old,new)
# Melt data table for easy operation
psl.ECE=melt(psl.ECE,id=c("lat","lon","startdate","targetdate"),variable.name="member",value.name="apsl")  

psl.ECE$member=as.numeric(as.character(psl.ECE$member))

#________________________________________
# Perform seasonal averages              \______________________________________________

sst.ECE=sst.ECE[,s.asst := ave(asst),by=.(startdate,member,lat,lon)]
# Eliminate monthly data
sst.ECE$asst = NULL
setnames(sst.ECE,"s.asst","asst")
# Remove repeated rows
sst.ECE = unique(sst.ECE, by=c("lat","lon","startdate","member"))

psl.ECE=psl.ECE[,s.apsl := ave(apsl),by=.(startdate,member,lat,lon)]
# Eliminate monthly data
psl.ECE$apsl = NULL
setnames(psl.ECE,"s.apsl","apsl")
# Remove repeated rows
psl.ECE = unique(psl.ECE, by=c("lat","lon","startdate","member"))

# If removing trend, here I should detrend asst by (lat,lon,member) after ordering by targetdate/startdate (startdate & targetdate by now are always unique w/e/other because lead is selected)

if(remove.trend==TRUE){
  # Remove linear trend of SST anomalies
  sst.ECE=sst.ECE[order(targetdate),]
  sst.ECE=sst.ECE[,dt.asst := detrend(asst),by=.(member,lat,lon)]
  sst.ECE$asst=NULL
  setnames(sst.ECE,"dt.asst","asst")
  
  # Remove linear trend of PSL anomalies
  psl.ECE=psl.ECE[order(targetdate),]
  psl.ECE=psl.ECE[,dt.apsl := detrend(apsl),by=.(member,lat,lon)]
  psl.ECE$apsl=NULL
  setnames(psl.ECE,"dt.apsl","apsl")
}

# Compute lead-dependent seasonal average for the ensemble mean
# SST
s.clim.sst.ECE=sst.ECE[member==99,]
s.clim.sst.ECE=s.clim.sst.ECE[,.(lat,lon,startdate,targetdate,asst)]
s.clim.sst.ECE=s.clim.sst.ECE[,s.clim := ave(asst),by=.(lat,lon)] #startdate & targetdate by now are always unique w/e/other because lead is selected
s.clim.sst.ECE=s.clim.sst.ECE[,.(lat,lon,s.clim)]
s.clim.sst.ECE=unique(s.clim.sst.ECE)
# Merge and compute seasonal anomaly
sst.ECE = merge(sst.ECE,s.clim.sst.ECE)
# Compute anomaly
sst.ECE=sst.ECE[,sa.asst := (asst - s.clim)]
sst.ECE$asst=NULL
sst.ECE$s.clim=NULL
setnames(sst.ECE,"sa.asst","asst")
# SLP
s.clim.psl.ECE=psl.ECE[member==99,]
s.clim.psl.ECE=s.clim.psl.ECE[,.(lat,lon,startdate,targetdate,apsl)]
s.clim.psl.ECE=s.clim.psl.ECE[,s.clim := ave(apsl),by=.(lat,lon)] #startdate & targetdate by now are always unique w/e/other because lead is selected
s.clim.psl.ECE=s.clim.psl.ECE[,.(lat,lon,s.clim)]
s.clim.psl.ECE=unique(s.clim.psl.ECE)
# Merge and compute seasonal anomaly
psl.ECE = merge(psl.ECE,s.clim.psl.ECE)
# Compute anomaly
psl.ECE=psl.ECE[,sa.apsl := (apsl - s.clim)]
psl.ECE$apsl=NULL
psl.ECE$s.clim=NULL
setnames(psl.ECE,"sa.apsl","apsl")

rm("s.clim.psl.ECE","s.clim.sst.ECE")

#________________________________________
# Apply latitude weight to the anomalies \______________________________________________

sst.ECE$asst=sst.ECE$asst*sqrt(cos(sst.ECE$lat*pi/180))
psl.ECE$apsl=psl.ECE$apsl*sqrt(cos(psl.ECE$lat*pi/180))

#______________________________________________________________
# Rearrange a data table to use all 15 members to compute SVD  \_________________________

# First, eliminate ensemble mean
sst.ECE = sst.ECE[member<98,]
psl.ECE = psl.ECE[member<98,]

sst.ECE=sst.ECE[order(member),]
psl.ECE=psl.ECE[order(member),]

# Generate artificial different dates changing century according to member number
sst.ECE.allin1 = sst.ECE[,.(lat,lon,targetdate,member,asst)]
sst.ECE.allin1$date = sst.ECE.allin1$targetdate %m+% years(sst.ECE.allin1$member*100)
sst.ECE.allin1$targetdate=NULL
sst.ECE.allin1$member=NULL

psl.ECE.allin1 = psl.ECE[,.(lat,lon,targetdate,member,apsl)]
psl.ECE.allin1$date = psl.ECE.allin1$targetdate %m+% years(psl.ECE.allin1$member*100)
psl.ECE.allin1$targetdate=NULL
psl.ECE.allin1$member=NULL

#________________________________________
# Remove NA values                       \______________________________________________

sst.ECE.allin1 = sst.ECE.allin1[!is.na(asst),]
psl.ECE.allin1 = psl.ECE.allin1[!is.na(apsl),]

##############################################
#
#   SVD - Atlantic SST & SLP
#
##############################################

# Change data to rows=space, columns=time
S = makematrix(sst.ECE.allin1,"asst ~ lat + lon | date")
P = makematrix(psl.ECE.allin1,"apsl ~ lat + lon | date")
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
exp.coef = data.table(date = unique(sst.ECE.allin1$date),ec1.sst=A[,1],ec2.sst=A[,2],ec3.sst=A[,3],ec1.slp=B[,1],ec2.slp=B[,2],ec3.slp=B[,3])

#***********************
# Rearrange data and separate by member, using the date transformation inversed
lastdate=max(sst.ECE$targetdate)

exp.coef$member=1
# exp.coef[date>(lastdate %m+% years(1*100)),]$member=2
for(i in 1:15){
  exp.coef[eval(parse(text=paste0("date>(as.Date(lastdate) %m+% years(",(i-1),"*100))"))),]$member=(i)
}
exp.coef$realdate=exp.coef$date %m-% years(exp.coef$member*100)
exp.coef$date=NULL
setnames(exp.coef,"realdate","targetdate")


#------------------------
# SST
exp.coef.sst1=exp.coef[,.(targetdate,member,ec1.sst)]

# Normalize ECs according to its sd
exp.coef.sst1=exp.coef.sst1[,ec1.sst.norm := ec1.sst/sd(ec1.sst,na.rm=TRUE),by="member"] # changed accordingly not to use smoothed version
# Merge with asst data to compute homogeneous correlation maps:
sst.ECE.allin1 = merge(sst.ECE,exp.coef.sst1[,.(targetdate,member,ec1.sst,ec1.sst.norm)],by=c("targetdate","member"))
sst.ECE.allin1 = sst.ECE.allin1[,hcm.1 := cor(asst,ec1.sst),by=.(lat,lon)] #All members together CORRELATION WITH SMOOTHED EC OR WITH RAW EC?

#------------------------
# SLP
exp.coef.slp1=exp.coef[,.(targetdate,member,ec1.slp)]

# Normalize ECs according to its sd
exp.coef.slp1=exp.coef.slp1[,ec1.slp.norm := ec1.slp/sd(ec1.slp,na.rm=TRUE),by="member"]  # changed accordingly not to use smoothed version

# Merge with aslp data to compute homogeneous correlation maps:
psl.ECE.allin1 = merge(psl.ECE,exp.coef.slp1[,.(targetdate,member,ec1.slp,ec1.slp.norm)],by=c("targetdate","member"))
psl.ECE.allin1 = psl.ECE.allin1[,hcm.1 := cor(apsl,ec1.slp),by=.(lat,lon)] #All members together CORRELATION WITH SMOOTHED EC OR WITH RAW EC?

#------------------------
rm(A,B,C,SV,SV2,S,P)

# Merge both expansion coefficients of SVD mode#1
exp.coef.separated=merge(exp.coef.sst1,exp.coef.slp1)

# Compute max,p80,median,p20,min by member
exp.coef.separated[,sst.min := min(ec1.sst.norm),by="targetdate"]
exp.coef.separated[,sst.p20 := quantile(ec1.sst.norm,0.13333,na.rm=TRUE),by="targetdate"]
exp.coef.separated[,sst.med := mean(ec1.sst.norm),by="targetdate"]
exp.coef.separated[,sst.p80 := quantile(ec1.sst.norm,0.86667,na.rm=TRUE),by="targetdate"]
exp.coef.separated[,sst.max := max(ec1.sst.norm),by="targetdate"]

exp.coef.separated[,slp.min := min(ec1.slp.norm),by="targetdate"]
exp.coef.separated[,slp.p20 := quantile(ec1.slp.norm,0.13333,na.rm=TRUE),by="targetdate"]
exp.coef.separated[,slp.med := mean(ec1.slp.norm),by="targetdate"]
exp.coef.separated[,slp.p80 := quantile(ec1.slp.norm,0.86667,na.rm=TRUE),by="targetdate"]
exp.coef.separated[,slp.max := max(ec1.slp.norm),by="targetdate"]  

#________________________________________
# Plotting and saving                    \______________________________________________

# Change longitude to -180:180 to plot correctly
psl.ECE.allin1[lon>180]$lon=psl.ECE.allin1[lon>180]$lon-360
sst.ECE.allin1[lon>180]$lon=sst.ECE.allin1[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-0.9
bmax=0.9
# bstep=0.1
bbreaks.contours=c(-99,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,99)
breaks.dates=exp.coef.separated[member==1,]$targetdate
limits.EC=2.5

# Need to select a single date as values are repeated by date
g1 = ggplot() +
  geom_contour_fill(data=sst.ECE.allin1[targetdate==sst.ECE.allin1$targetdate[1] & member==1],aes(lon, lat, z = hcm.1, fill=stat(level)),breaks=bbreaks.contours,na.fill=TRUE)+
  # scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
  #                      breaks=bbreaks.cbar,
  #                      limits=c(bmin,bmax),
  #                      guide = guide_colorstrip(),
  #                      labels=labels.cbar,
  #                      oob  = scales::squish)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       limits=c(bmin,bmax),
                       super = ScaleDiscretised,
                       guide = guide_colorsteps(),
                       oob  = scales::squish)+
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 10))+
  new_scale_color() +
  geom_contour(data=psl.ECE.allin1[targetdate==psl.ECE.allin1$targetdate[1] & member==1],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=psl.ECE.allin1[targetdate==psl.ECE.allin1$targetdate[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("SVD1 as homogeneous correlation map: SST (shaded), SLP (contours) SCF: ",as.character(round(SCF[1],1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.separated[member==1,],aes(targetdate, sst.med),col="#de77ae",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated[member==1,], aes(targetdate, ymin=sst.p20 , ymax=sst.p80 ),fill="#de77ae",alpha=0.4)+
  # geom_line(data=exp.coef.separated,aes(date, sst.max),col="red",alpha=0.4,size=0.2)+
  # geom_line(data=exp.coef.separated,aes(date, sst.min),col="red",alpha=0.4,size=0.2)+
  
  geom_line(data=exp.coef.separated[member==1,],aes(targetdate, slp.med),col="#7fbc41",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated[member==1,], aes(targetdate, ymin=slp.p20 , ymax=slp.p80 ),fill="#7fbc41",alpha=0.3)+
  # geom_line(data=exp.coef.separated,aes(date, slp.max),col="#7fbc41",alpha=0.4,size=0.2)+
  # geom_line(data=exp.coef.separated,aes(date, slp.min),col="#7fbc41",alpha=0.4,size=0.2)  
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),4)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(ceiling(-limits.EC),floor(limits.EC),1),limits=c(-limits.EC,limits.EC),expand = c(0, 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("targetdate")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("EC1 SST (pink), SLP (green). r=",as.character(round(cor(exp.coef$ec1.sst,exp.coef$ec1.slp,use = "pairwise.complete.obs"),2))," for raw EC. Shading: [P13,P86] interval (9 less extreme mmbs)"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))

# Density plot of expansion coefficients
g3 <- ggplot(data=exp.coef.separated, aes(x=ec1.sst.norm, y=ec1.slp.norm)) +
  theme_bw()+
  geom_hex(bins=70) +
  scale_fill_continuous(type="viridis") +
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 6))+
  scale_x_continuous(limits=c(-4.5,4.5),expand = c(0., 0.))+
  scale_y_continuous(limits=c(-4.5,4.5),expand = c(0., 0.))+
  geom_abline(intercept = 0,color="red")+
  xlab("EC1: SST (norm. units)")+ylab("EC1: SLP (norm. units)")+
  ggtitle(paste0("Density plot of norm. EC (r=",as.character(round(cor(exp.coef.separated$ec1.sst.norm,exp.coef.separated$ec1.slp.norm,use = "pairwise.complete.obs"),2))," )"))+
  theme(axis.text = element_text(size=11),axis.title = element_text(size=11),title = element_text(size=10))

# Save

if(remove.trend==TRUE){
    fig <- grid.arrange(g1,g2,g3, layout_matrix=rbind(c(1,1,4),c(2,2,3)),top = textGrob(paste0(sel.season,", DCPP, lead year: ",sel.lead,"; SVD of SST-SLP anomalies (no trend, weighted) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_dcpp_SVD_SST_SLP_weighted_notrend_allmembers_lead_year_",sel.lead,".png"),plot=fig,width = 12, height = 8)
    
    # Save expansion coefficients & patterns
    save(sst.ECE.allin1,psl.ECE.allin1,exp.coef.separated,file=paste0("/home/maralv/data/ECs_SVD1_",sel.season,"_dcpp_weighted_notrend_allmembers_lead_year_",sel.lead,".rda"))
    

}else{ # No trend

    fig <- grid.arrange(g1,g2,g3, layout_matrix=rbind(c(1,1,4),c(2,2,3)),top = textGrob(paste0(sel.season,", DCPP, lead year: ",sel.lead,"; SVD of SST-SLP anomalies ( weighted) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/",sel.season,"_dcpp_SVD_SST_SLP_weighted_allmembers_lead_year_",sel.lead,".png"),plot=fig,width = 12, height = 8)

    # Save expansion coefficients & patterns
    save(sst.ECE.allin1,psl.ECE.allin1,exp.coef.separated,file=paste0("/home/maralv/data/ECs_SVD1_",sel.season,"_dcpp_weighted_allmembers_lead_year_",sel.lead,".rda"))
} #endif trend

rm(list=ls())
graphics.off()

} #end for leads