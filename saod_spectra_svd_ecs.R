# This script computes the spectra from expansion coefficients
#
# M. Alvarez (2020)
#--------------------------------------------------------------------------------
rm(list=ls())
graphics.off()

setwd("/home/maralv/")

library("metR")
library('data.table')
library("ggplot2")
library("gridExtra")
library("grid")

#---------------------------------------------------------------------------------------
#  functions 
#---------------------------------------------------------------------------------------

library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#---------------------------------------------------------------------------------------
#  Functions to compute power spectra by Elio Campitelli (GitHub @eliocamp)
#---------------------------------------------------------------------------------------
#' FFT spectrum
#'
#' A light wrapper around [stats::spec.pgram()] and [stats::ar()]. Returns the spectrum of
#' the signal and the spectrum of the fitted autoregressiv model.
#'
#' @param x numeric vector
#' @param B number of bootstrap samples.
#' @param spans vector of odd integers giving the widths of modified Daniell
#' smoothers to be used to smooth the periodogram
#' @param ... other arguments passed to [stats::spec.pgram()]
#'
#' @export
fftspectrum <- function(x, spans = NULL, B = 10000, ..., probs = 0.95) {
  mtm <- spec.pgram(x, spans = spans, ..., plot = FALSE)  
  out <- as.data.table(mtm[c("freq", "spec")])  
  ar <- ar(ts(x))
  # rho <- a$ar
  # var <- a$var.pred
  out[, ar_spectrum := arspectrum(mtm$freq, ar$ar, ar$var.pred)]
  out[, c(scales::percent(probs)) := null_ar_spectrum(B = B, length(x), ar, spans = spans, ..., probs = probs)]  
  return(out[])
}

arspectrum <- function(freq, rho, var) {
  k <- seq_along(rho)
  e <- vapply(freq, function(f) sum(rho * exp(-2*1i*pi*f*k)), complex(1))
  var / (Mod(1 - e))^2
}

null_ar_spectrum_ <- function(B = 100, n, ar, spans = NULL, ..., probs = 0.95) {
  y <- as.vector(arima.sim(model = list(ar = ar$ar), n = n))*sqrt(ar$var.pred)
  nfreq <- length(spec.pgram(y, spans = spans, ..., plot = FALSE)$spec)
  boots <- vapply(seq_len(B), function(b) {
    y <- as.vector(arima.sim(model = list(ar = ar$ar), n = n))*sqrt(ar$var.pred)
    spec.pgram(y, spans = spans, plot = FALSE)$spec  
    
  }, numeric(nfreq))
  data.table::as.data.table((apply(boots, 1, quantile, probs = probs)))
}

null_ar_spectrum <- memoise::memoise(null_ar_spectrum_)
# null_ar_spectrum <- null_ar_spectrum_

null_spec <- memoise::memoise(function(x, spans, B = 1000, ..., probs = 0.95) {
  b <- boot::boot(x, function(d, i)  spec.pgram(d[i], spans = spans,
                                                ...,
                                                plot = FALSE)$spec,
                  R = B)  
  apply(b$t, 2, quantile, probs = probs)
})

#---------------------------------------------------------------------------------------
#  Man Program
#---------------------------------------------------------------------------------------

# Load Expansion Coefficients

load(file="SAO_monthly_expansion_coefs.RData")

#-------
# SVD1

out = fftspectrum(exp.coef$ec.sst.1, spans = c(3,5), B = 10000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

# Plot

g1 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#238443",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#78c679",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#d9f0a3",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SAOD (SVD1) - SST Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

out = fftspectrum(exp.coef$ec.slp.1, spans = c(3,5), B = 10000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

# Plot

g2 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#ae017e",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#f768a1",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#fa9fb5",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SAOD (SVD1) - SLP Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

#-------
# SVD2

out = fftspectrum(exp.coef$ec.sst.2, spans = c(3,5), B = 10000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

# Plot

g3 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#238443",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#78c679",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#d9f0a3",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SVD2 - SST Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

out = fftspectrum(exp.coef$ec.slp.2, spans = c(3,5), B = 10000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

# Plot

g4 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#ae017e",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#f768a1",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#fa9fb5",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SVD2 - SLP Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

#-------
# SVD3

out = fftspectrum(exp.coef$ec.sst.3, spans = c(3,5), B = 10000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

# Plot

g5 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#238443",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#78c679",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#d9f0a3",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SVD3 - SST Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

out = fftspectrum(exp.coef$ec.slp.3, spans = c(3,5), B = 10000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

# Plot

g6 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#ae017e",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#f768a1",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#fa9fb5",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SVD3 - SLP Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

fig <- grid.arrange(g1,g2,g3,g4,g5,g6, ncol = 2,top = textGrob(paste0("SST-SLP monthly SVD Expansion Coefficients Spectra"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/spectra.png"),plot=fig,width = 13, height = 12)


#######################3
# Performs spectra of smooth (but not normalized) coefficients

#########################################################################################
# SST

library('data.table')

data1 = as.data.table(exp.coef$ec.sst.1)
data1 = frollmean(data1,13,align="center",fill=NA)[[1]]
# remove NA
data1=data1[!is.na(data1)]

out = fftspectrum(data1, spans = c(3,5), B = 1000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

# Plot

g7 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#ae017e",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#f768a1",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#fa9fb5",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SAOD (SVD1) - SST Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

data1 = as.data.table(exp.coef$ec.slp.1)
data1 = frollmean(data1,13,align="center",fill=NA)[[1]]
# remove NA
data1=data1[!is.na(data1)]

out = fftspectrum(data1, spans = c(3,5), B = 1000, probs = 0.90)

out$period = 1/out$freq
out$periodYY = out$period/12

g8 <- ggplot()+
  theme_bw()+
  geom_line(data=out,aes(periodYY,spec),col="#238443",alpha=0.9)+
  geom_line(data=out,aes(periodYY,ar_spectrum),col="#78c679",alpha=0.9,linetype = "dashed")+
  geom_ribbon(data=out, aes(periodYY, ymin=(spec*0) , ymax=`90%` ),fill="#d9f0a3",alpha=0.5)+
  scale_x_continuous(trans=reverselog_trans(10),breaks=c(72,30,20,14,9,5,3,2,1))+
  labs(x="Period (years)",y="Power",title = "SAOD (SVD1) - SLP Exp. Coeff. Spectrum (90%)")+
  theme(text = element_text(size=14),title = element_text(size=14),axis.text = element_text(size = 14))

fig <- grid.arrange(g7,g8, ncol = 2)
ggsave(filename=paste0("/home/maralv/Dropbox/DMI/Figures/spectra_SVD1_ECs_observations_monthly.png"),plot=fig,width = 13, height = 6)
