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
  
  #________________________________________
  # Linear trend removal while monthly     \______________________________________________
  
  # For now, only for ensemble mean respect to itself
  
  # SST
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
  
  # SLP
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
  # Apply latitude weight to the anomalies \______________________________________________
  
  # SST
  for (mmb in members) {
    eval(parse(text = paste0("sst.ECE$asst.",mmb, "=","sst.ECE$asst.",mmb,"*sqrt(cos(sst.ECE$lat*pi/180))" )))
  }
  sst.ECE$asst.em=sst.ECE$asst.em*sqrt(cos(sst.ECE$lat*pi/180))
  
  # SLP
  for (mmb in members) {
    eval(parse(text = paste0("psl.ECE$apsl.",mmb, "=","psl.ECE$apsl.",mmb,"*sqrt(cos(psl.ECE$lat*pi/180))" )))
  }
  psl.ECE$apsl.em=psl.ECE$apsl.em*sqrt(cos(psl.ECE$lat*pi/180))
  
  
  #______________________________________________________________
  # Rearrange a data table to use all 15 members to compute EOF  \_________________________
  
  #### SST
  # Create a lat-lon-targetdase base data table
  sst.ECE.allin1 = sst.ECE[,c("lat","lon","targetdate")]
  sst.ECE.2add = sst.ECE.allin1
  
  # replicate rows 14 times (as there are 15 members to be used) 
  for(reps in 1:14){
    sst.ECE.2add$targetdate=sst.ECE.2add$targetdate+(100*365) # Add 100 years to each round
    sst.ECE.allin1 = rbind(sst.ECE.allin1,sst.ECE.2add)
  }
  
  sst.ECE.allin1$asst=NA_real_
  
  l=length(sst.ECE$targetdate)
  
  sst.ECE.allin1$asst[1:l]=sst.ECE$asst.1
  sst.ECE.allin1$asst[(1+l*1):(l*2)]=sst.ECE$asst.2
  sst.ECE.allin1$asst[(1+l*2):(l*3)]=sst.ECE$asst.4
  sst.ECE.allin1$asst[(1+l*3):(l*4)]=sst.ECE$asst.5
  sst.ECE.allin1$asst[(1+l*4):(l*5)]=sst.ECE$asst.6
  sst.ECE.allin1$asst[(1+l*5):(l*6)]=sst.ECE$asst.7
  sst.ECE.allin1$asst[(1+l*6):(l*7)]=sst.ECE$asst.8
  sst.ECE.allin1$asst[(1+l*7):(l*8)]=sst.ECE$asst.9
  sst.ECE.allin1$asst[(1+l*8):(l*9)]=sst.ECE$asst.10
  sst.ECE.allin1$asst[(1+l*9):(l*10)]=sst.ECE$asst.11
  sst.ECE.allin1$asst[(1+l*10):(l*11)]=sst.ECE$asst.12
  sst.ECE.allin1$asst[(1+l*11):(l*12)]=sst.ECE$asst.13
  sst.ECE.allin1$asst[(1+l*12):(l*13)]=sst.ECE$asst.14
  sst.ECE.allin1$asst[(1+l*13):(l*14)]=sst.ECE$asst.15
  sst.ECE.allin1$asst[(1+l*14):(l*15)]=sst.ECE$asst.16
  
  
  #### PSL
  # Create a lat-lon-targetdase base data table
  psl.ECE.allin1 = psl.ECE[,c("lat","lon","targetdate")]
  psl.ECE.2add = psl.ECE.allin1
  
  # replicate rows 14 times (as there are 15 members to be used) 
  for(reps in 1:14){
    psl.ECE.2add$targetdate=psl.ECE.2add$targetdate+(100*365) # Add 100 years to each round
    psl.ECE.allin1 = rbind(psl.ECE.allin1,psl.ECE.2add)
  }
  
  psl.ECE.allin1$apsl=NA_real_
  
  l=length(psl.ECE$targetdate)
  
  psl.ECE.allin1$apsl[1:l]=psl.ECE$apsl.1
  psl.ECE.allin1$apsl[(1+l*1):(l*2)]=psl.ECE$apsl.2
  psl.ECE.allin1$apsl[(1+l*2):(l*3)]=psl.ECE$apsl.4
  psl.ECE.allin1$apsl[(1+l*3):(l*4)]=psl.ECE$apsl.5
  psl.ECE.allin1$apsl[(1+l*4):(l*5)]=psl.ECE$apsl.6
  psl.ECE.allin1$apsl[(1+l*5):(l*6)]=psl.ECE$apsl.7
  psl.ECE.allin1$apsl[(1+l*6):(l*7)]=psl.ECE$apsl.8
  psl.ECE.allin1$apsl[(1+l*7):(l*8)]=psl.ECE$apsl.9
  psl.ECE.allin1$apsl[(1+l*8):(l*9)]=psl.ECE$apsl.10
  psl.ECE.allin1$apsl[(1+l*9):(l*10)]=psl.ECE$apsl.11
  psl.ECE.allin1$apsl[(1+l*10):(l*11)]=psl.ECE$apsl.12
  psl.ECE.allin1$apsl[(1+l*11):(l*12)]=psl.ECE$apsl.13
  psl.ECE.allin1$apsl[(1+l*12):(l*13)]=psl.ECE$apsl.14
  psl.ECE.allin1$apsl[(1+l*13):(l*14)]=psl.ECE$apsl.15
  psl.ECE.allin1$apsl[(1+l*14):(l*15)]=psl.ECE$apsl.16
  
  
  rm("sst.ECE.2add","psl.ECE.2add")
  
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
  
  # Rename targetdate
  setnames(sst.ECE.allin1,"targetdate","date")
  setnames(psl.ECE.allin1,"targetdate","date")
  
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
  exp.coef = data.table(date = unique(sst.ECE.allin1$date),ec.sst.1=A[,1],ec.sst.2=A[,2],ec.sst.3=A[,3],ec.slp.1=B[,1],ec.slp.2=B[,2],ec.slp.3=B[,3])
  
  # Create data table with smoothed expansion coefficients and normalize respect to sd
  exp.coef.norm = data.table(date = unique(sst.ECE.allin1$date),ec.sst.1=rollmean(A[,1],13,align="center",fill=NA),ec.sst.2=rollmean(A[,2],13,align="center",fill=NA),ec.sst.3=rollmean(A[,3],13,align="center",fill=NA),ec.slp.1=rollmean(B[,1],13,align="center",fill=NA),ec.slp.2=rollmean(B[,2],13,align="center",fill=NA),ec.slp.3=rollmean(B[,3],13,align="center",fill=NA))
  #exp.coef.norm = data.table(date = unique(sst.ECE$date),ec.sst.1=A[,1],ec.sst.2=A[,2],ec.sst.3=A[,3],ec.slp.1=B[,1],ec.slp.2=B[,2],ec.slp.3=B[,3])
  exp.coef.norm = exp.coef.norm[,.(date,ec.sst.1=ec.sst.1/sd(ec.sst.1,na.rm=TRUE),ec.sst.2=ec.sst.2/sd(ec.sst.2,na.rm=TRUE),ec.sst.3=ec.sst.3/sd(ec.sst.3,na.rm=TRUE),ec.slp.1=ec.slp.1/sd(ec.slp.1,na.rm=TRUE),ec.slp.2=ec.slp.2/sd(ec.slp.2,na.rm=TRUE),ec.slp.3=ec.slp.3/sd(ec.slp.3,na.rm=TRUE))]
  
  rm(A,B,C,SV,SV2,S,P)
  
  # As I will present spatial patterns as homogeneous correlation maps:
  psl.ECE.allin1 = merge(psl.ECE.allin1,exp.coef[,.(date,ec.slp.1,ec.slp.2,ec.slp.3)],by="date")
  sst.ECE.allin1 = merge(sst.ECE.allin1,exp.coef[,.(date,ec.sst.1,ec.sst.2,ec.sst.3)],by="date")
  
  psl.ECE.allin1 = psl.ECE.allin1[,.(date,apsl,ec.slp.1,ec.slp.2,ec.slp.3,hcm.1 = cor(apsl,ec.slp.1),hcm.2 = cor(apsl,ec.slp.2),hcm.3 = cor(apsl,ec.slp.3)),by=.(lat,lon)]
  sst.ECE.allin1 = sst.ECE.allin1[,.(date,asst,ec.sst.1,ec.sst.2,ec.sst.3,hcm.1 = cor(asst,ec.sst.1),hcm.2 = cor(asst,ec.sst.2),hcm.3 = cor(asst,ec.sst.3)),by=.(lat,lon)]
  
  # Separate expansion coefficients
  # This should be done prior to applying moving average for smoothing
  my=max((sst.ECE$targetdate))
  a=sum(((exp.coef.norm$date)==my)*seq(1,length(exp.coef.norm$date),1))
  
  # Arrange in a data table to melt and then lead to easier operation
  # SST 1
  exp.coef.separated.sst1=exp.coef.norm[1:a,.(date,ec.sst.1)]
  setnames(exp.coef.separated.sst1,"ec.sst.1","1")
  exp.coef.separated.sst1[,`:=` ("2"=exp.coef.norm$ec.sst.1[(1+a*1):(a*2)])]
  for(i in c(4,5,6,7,8,9,10,11,12,13,14,15,16)){
    exp.coef.separated.sst1[,eval(parse(text=paste0("`:=` ('",i,"'=exp.coef.norm$ec.sst.1[(1+a*",(i-2),"):(a*",(i-1),")])"))) ]
  }
  exp.coef.separated.sst1=melt(exp.coef.separated.sst1,id="date",variable.name="member",value.name="ec1.sst")
  
  # SLP 1
  exp.coef.separated.slp1=exp.coef.norm[1:a,.(date,ec.slp.1)]
  setnames(exp.coef.separated.slp1,"ec.slp.1","1")
  exp.coef.separated.slp1[,`:=` ("2"=exp.coef.norm$ec.slp.1[(1+a*1):(a*2)])]
  for(i in c(4,5,6,7,8,9,10,11,12,13,14,15,16)){
    exp.coef.separated.slp1[,eval(parse(text=paste0("`:=` ('",i,"'=exp.coef.norm$ec.slp.1[(1+a*",(i-2),"):(a*",(i-1),")])"))) ]
  }
  exp.coef.separated.slp1=melt(exp.coef.separated.slp1,id="date",variable.name="member",value.name="ec1.slp")
  
  # Merge both expansion coefficients of SVD mode#1
  exp.coef.separated=merge(exp.coef.separated.sst1,exp.coef.separated.slp1)
  
  # Compute max,p80,median,p20,min by member
  exp.coef.separated[,sst.min := min(ec1.sst),by="date"]
  exp.coef.separated[,sst.p20 := quantile(ec1.sst,0.13333,na.rm=TRUE),by="date"]
  exp.coef.separated[,sst.med := mean(ec1.sst),by="date"]
  exp.coef.separated[,sst.p80 := quantile(ec1.sst,0.86667,na.rm=TRUE),by="date"]
  exp.coef.separated[,sst.max := max(ec1.sst),by="date"]
  
  exp.coef.separated[,slp.min := min(ec1.slp),by="date"]
  exp.coef.separated[,slp.p20 := quantile(ec1.slp,0.13333,na.rm=TRUE),by="date"]
  exp.coef.separated[,slp.med := mean(ec1.slp),by="date"]
  exp.coef.separated[,slp.p80 := quantile(ec1.slp,0.86667,na.rm=TRUE),by="date"]
  exp.coef.separated[,slp.max := max(ec1.slp),by="date"]  
  
  # setnames(exp.coef.separated.sst,"ec.sst.2","2.1")
  # setnames(exp.coef.separated.sst,"ec.sst.3","3.1")
  # 
  # aux=melt(exp.coef.separated.sst,id="date")
  # 
  # exp.coef.separated.sst[,`:=` (ec.sst.1.2=exp.coef.norm$ec.sst.1[(1+a*1):(a*2)], ec.sst.2.2=exp.coef.norm$ec.sst.2[(1+a*1):(a*2)], ec.sst.3.2=exp.coef.norm$ec.sst.3[(1+a*1):(a*2)])]
  # 
  # 
  # exp.coef.separated.slp=exp.coef.norm[1:a,]
  # setnames(exp.coef.separated.slp,"ec.slp.1","ec.slp.1.1")
  # setnames(exp.coef.separated.slp,"ec.slp.2","ec.slp.2.1")
  # setnames(exp.coef.separated.slp,"ec.slp.3","ec.slp.3.1")
  # 
  # exp.coef.separated[,`:=` (ec.slp.1.2=exp.coef.norm$ec.slp.1[(1+a*1):(a*2)], ec.slp.2.2=exp.coef.norm$ec.slp.2[(1+a*1):(a*2)], ec.slp.3.2=exp.coef.norm$ec.slp.3[(1+a*1):(a*2)])]
  #
    # for(i in c(4,5,6,7,8,9,10,11,12,13,14,15,16)){
  #   exp.coef.separated[,eval(parse(text=paste0("`:=` (ec.sst.1.",i,"=exp.coef.norm$ec.sst.1[(1+a*",(i-2),"):(a*",(i-1),")],  ec.sst.2.",i,"=exp.coef.norm$ec.sst.2[(1+a*",(i-2),"):(a*",(i-1),")],  ec.sst.3.",i,"=exp.coef.norm$ec.sst.3[(1+a*",(i-2),"):(a*",(i-1),")])"))) ]
  #   exp.coef.separated[,eval(parse(text=paste0("`:=` (ec.slp.1.",i,"=exp.coef.norm$ec.slp.1[(1+a*",(i-2),"):(a*",(i-1),")],  ec.slp.2.",i,"=exp.coef.norm$ec.slp.2[(1+a*",(i-2),"):(a*",(i-1),")],  ec.slp.3.",i,"=exp.coef.norm$ec.slp.3[(1+a*",(i-2),"):(a*",(i-1),")])"))) ]
  #   
  # }
  



#________________________________________
# Plotting and saving                    \______________________________________________

# Change longitude to -180:180 to plot correctly
psl.ECE.allin1[lon>180]$lon=psl.ECE.allin1[lon>180]$lon-360
sst.ECE.allin1[lon>180]$lon=sst.ECE.allin1[lon>180]$lon-360
map.world <- map_data ("world2", wrap = c(-180,180))

bmin=-0.9
bmax=0.9
bstep=0.1
bbreaks.contours=c(-99,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,99)
bbreaks.cbar=c(-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
labels.cbar=as.character(bbreaks.cbar)
labels.cbar[1]=""
labels.cbar[length(labels.cbar)]=""
breaks.dates=exp.coef.separated[member==1,]$date


# Need to select a single date as values are repeated by date
g1 = ggplot() +
  geom_contour_fill(data=sst.ECE.allin1[date==sst.ECE.allin1$date[1]],aes(lon, lat, z = hcm.1),breaks=bbreaks.contours,na.fill=TRUE)+
  scale_fill_distiller(name="SST",palette="RdBu",direction=-1,
                       breaks=bbreaks.cbar,
                       limits=c(bmin,bmax),
                       guide = guide_colorstrip(),
                       labels=labels.cbar,
                       oob  = scales::squish)+
  new_scale_color() +
  geom_contour(data=psl.ECE.allin1[date==psl.ECE.allin1$date[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),color="black",size=0.25)+
  geom_text_contour(data=psl.ECE.allin1[date==psl.ECE.allin1$date[1]],aes(lon, lat, z = hcm.1),breaks=seq(-1,1,0.1),stroke = 0.1,min.size = 10)+
  scale_x_longitude(breaks=seq(-70,20,20))+
  scale_y_latitude(breaks=seq(-40,0,10))+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black", inherit.aes = F)+
  ggtitle(paste0("SVD1 as homogeneous correlation map: SST (shaded), SLP (contours) SCF: ",as.character(round(SCF[1],1)),"%"))+
  theme(axis.text=element_text(size=12),title = element_text(size=10))

g2 = ggplot() +
  theme_bw()+
  geom_line(data=exp.coef.separated,aes(date, sst.med),col="#de77ae",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated, aes(date, ymin=sst.p20 , ymax=sst.p80 ),fill="#de77ae",alpha=0.4)+
  # geom_line(data=exp.coef.separated,aes(date, sst.max),col="red",alpha=0.4,size=0.2)+
  # geom_line(data=exp.coef.separated,aes(date, sst.min),col="red",alpha=0.4,size=0.2)+
  
  geom_line(data=exp.coef.separated,aes(date, slp.med),col="#7fbc41",alpha=0.8,size=0.7)+
  geom_ribbon(data=exp.coef.separated, aes(date, ymin=slp.p20 , ymax=slp.p80 ),fill="#7fbc41",alpha=0.3)+
  # geom_line(data=exp.coef.separated,aes(date, slp.max),col="#7fbc41",alpha=0.4,size=0.2)+
  # geom_line(data=exp.coef.separated,aes(date, slp.min),col="#7fbc41",alpha=0.4,size=0.2)  
  
  scale_x_date(breaks=breaks.dates[seq(1,length(breaks.dates),48)],date_labels = "%Y",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1),limits=c(-3,3),expand = c(0., 0.))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Date")+ylab("EC1 (normalized units)")+
  ggtitle(paste0("Normalized Expansion Coefficients: SST (pink), SLP (green). r=",as.character(round(cor(exp.coef$ec.sst.1,exp.coef$ec.slp.1,use = "pairwise.complete.obs"),2))," for raw EC. Shading denotes [P13,P86] interval (9 less extreme members)"))+
  theme(axis.title = element_text(size=11),title = element_text(size=10))




# Save Figure and PCs for the lead

if(remove.trend==TRUE){

    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0("All months , historical: SVD of SST-SLP anomalies (no trend, weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    # ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/AllMonths_historical_SVD_SST_SLP_weighted_notrend_ensmean.png"),plot=fig,width = 8, height = 8)
    
  
}else{
  # save(sst.pcs,sst.eof,file=paste0("/home/maralv/data/AllMonths_historical_ECEarth3_PCs_EOFs_SST_weighted_allmembers.RData"))

    fig <- grid.arrange(g1,g2, ncol = 1,top = textGrob(paste0("All months, historical: SVD of SST-SLP anomalies (weighted by cos(lat)) EC-Earth3"),gp=gpar(fontsize=13,font=3)))
    # ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/AllMonths_historical_SVD_SST_SLP_weighted_ensmean.png"),plot=fig,width = 8, height = 8)

} 

