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
  ######################################################
  
  # Rearrange data using number of member as variable
  # STT
  # rename according to member
  old = c("asst.1","asst.2","asst.4","asst.5","asst.6","asst.7","asst.8","asst.9","asst.10","asst.11","asst.12","asst.13","asst.14","asst.15","asst.16","asst.em")
  new = c("1","2","4","5","6","7","8","9","10","11","12","13","14","15","16","99")
  setnames(sst.ECE,old,new)
  # Melt data table for easy operation
  sst.ECE=melt(sst.ECE,id=c("lat","lon","targetdate","targetmonth"),variable.name="member",value.name="asst")
  
  sst.ECE$member=as.numeric(sst.ECE$member)
  
  # SLP
  # rename according to member
  old = c("apsl.1","apsl.2","apsl.4","apsl.5","apsl.6","apsl.7","apsl.8","apsl.9","apsl.10","apsl.11","apsl.12","apsl.13","apsl.14","apsl.15","apsl.16","apsl.em")
  setnames(psl.ECE,old,new)
  # Melt data table for easy operation
  psl.ECE=melt(psl.ECE,id=c("lat","lon","targetdate","targetmonth"),variable.name="member",value.name="apsl")  
  
  psl.ECE$member=as.numeric(psl.ECE$member)
  
  #________________________________________
  # Linear trend removal while monthly     \______________________________________________
  
  # SST
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
  exp.coef$member=1
  # Last real date is "2014-12-16"
  lastdate="2014-12-16"
  exp.coef[date>(as.Date(lastdate) %m+% years(1*100)),]$member=2
  for(i in c(2,4,5,6,7,8,9,10,11,12,13,14,15,16)){
    exp.coef[eval(parse(text=paste0("date>(as.Date(lastdate) %m+% years(",i,"*100))"))),]$member=(i+1)
  }
  exp.coef$realdate=exp.coef$date %m-% years(exp.coef$member*100)
  exp.coef$date=NULL
  setnames(exp.coef,"realdate","date")
  
  #------------------------
  # SST
  exp.coef.sst1=exp.coef[,.(date,member,ec1.sst)]

  # Generate smoothed expansion coefficients and then normalize them according to its sd
  exp.coef.sst1=exp.coef.sst1[,ec1.sst.smth := rollmean(ec1.sst,13,align="center",fill=NA),by="member"]
  exp.coef.sst1=exp.coef.sst1[,ec1.sst.norm := ec1.sst.smth/sd(ec1.sst.smth,na.rm=TRUE),by="member"]
  
  # Merge with asst data to compute homogeneous correlation maps:
  setnames(sst.ECE,"targetdate","date")
  sst.ECE.allin1 = merge(sst.ECE,exp.coef.sst1[,.(date,member,ec1.sst,ec1.sst.norm)],by=c("date","member"))
  sst.ECE.allin1 = sst.ECE.allin1[,hcm.1 := cor(asst,ec1.sst),by=.(lat,lon)] #All members together CORRELATION WITH SMOOTHED EC OR WITH RAW EC?

  #------------------------
  # SLP
  exp.coef.slp1=exp.coef[,.(date,member,ec1.slp)]
  
  # Generate smoothed expansion coefficients and then normalize them according to its sd
  exp.coef.slp1=exp.coef.slp1[,ec1.slp.smth := rollmean(ec1.slp,13,align="center",fill=NA),by="member"]
  exp.coef.slp1=exp.coef.slp1[,ec1.slp.norm := ec1.slp.smth/sd(ec1.slp.smth,na.rm=TRUE),by="member"]
  
  # Merge with aslp data to compute homogeneous correlation maps:
  setnames(psl.ECE,"targetdate","date")
  psl.ECE.allin1 = merge(psl.ECE,exp.coef.slp1[,.(date,member,ec1.slp,ec1.slp.norm)],by=c("date","member"))
  psl.ECE.allin1 = psl.ECE.allin1[,hcm.1 := cor(apsl,ec1.slp),by=.(lat,lon)] #All members together CORRELATION WITH SMOOTHED EC OR WITH RAW EC?
  #------------------------
  rm(A,B,C,SV,SV2,S,P)
  
  # Merge both expansion coefficients of SVD mode#1
  exp.coef.separated=merge(exp.coef.sst1,exp.coef.slp1)
  
  # Compute max,p80,median,p20,min by member
  exp.coef.separated[,sst.min := min(ec1.sst.norm),by="date"]
  exp.coef.separated[,sst.p20 := quantile(ec1.sst.norm,0.13333,na.rm=TRUE),by="date"]
  exp.coef.separated[,sst.med := mean(ec1.sst.norm),by="date"]
  exp.coef.separated[,sst.p80 := quantile(ec1.sst.norm,0.86667,na.rm=TRUE),by="date"]
  exp.coef.separated[,sst.max := max(ec1.sst.norm),by="date"]
  
  exp.coef.separated[,slp.min := min(ec1.slp.norm),by="date"]
  exp.coef.separated[,slp.p20 := quantile(ec1.slp.norm,0.13333,na.rm=TRUE),by="date"]
  exp.coef.separated[,slp.med := mean(ec1.slp.norm),by="date"]
  exp.coef.separated[,slp.p80 := quantile(ec1.slp.norm,0.86667,na.rm=TRUE),by="date"]
  exp.coef.separated[,slp.max := max(ec1.slp.norm),by="date"]  

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
breaks.dates=exp.coef.separated[member==1,]$date

# Need to select a single date as values are repeated by date
g1 = ggplot() +
  geom_contour_fill(data=sst.ECE.allin1[date==sst.ECE.allin1$date[1] & member==1],aes(lon, lat, z = hcm.1, fill=stat(level)),breaks=bbreaks.contours,na.fill=TRUE)+
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
  geom_contour(data=psl.ECE.allin1[date==psl.ECE.allin1$date[1] & member==1],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=psl.ECE.allin1[date==psl.ECE.allin1$date[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("SVD1 as homogeneous correlation map: SST (shaded), SLP (contours) SCF: ",as.character(round(SCF[1],1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.separated[member==1,],aes(date, sst.med),col="#de77ae",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated[member==1,], aes(date, ymin=sst.p20 , ymax=sst.p80 ),fill="#de77ae",alpha=0.4)+
  # geom_line(data=exp.coef.separated,aes(date, sst.max),col="red",alpha=0.4,size=0.2)+
  # geom_line(data=exp.coef.separated,aes(date, sst.min),col="red",alpha=0.4,size=0.2)+
  
  geom_line(data=exp.coef.separated[member==1,],aes(date, slp.med),col="#7fbc41",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated[member==1,], aes(date, ymin=slp.p20 , ymax=slp.p80 ),fill="#7fbc41",alpha=0.3)+
  # geom_line(data=exp.coef.separated,aes(date, slp.max),col="#7fbc41",alpha=0.4,size=0.2)+
  # geom_line(data=exp.coef.separated,aes(date, slp.min),col="#7fbc41",alpha=0.4,size=0.2)  
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),48)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3,3),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("EC1 (normalized units)")+
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

# Save Figure and PCs for the lead

if(remove.trend==TRUE){

    fig <- grid.arrange(g1,g2,g3, layout_matrix=rbind(c(1,1,4),c(2,2,3)),top = textGrob(paste0("All months , historical: SVD of SST-SLP anomalies (no trend, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/AllMonths_historical_SVD_SST_SLP_weighted_notrend_allmembers_v2.png"),plot=fig,width = 12, height = 8)
  
}else{
  # save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/AllMonths_historical_ECEarth3_PCs_EOFs_SST_weighted_allmembers.RData"))

    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0("All months, historical: SVD of SST-SLP anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    # ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/AllMonths_historical_SVD_SST_SLP_weighted_ensmean.png"),plot=fig,width = 8, height = 8)

} 



 
