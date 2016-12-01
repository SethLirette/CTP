library(fields)
library(oro.dicom)
library(oro.nifti)
library(MASS)
library(zoo)
library(caTools)



#######################################################
### Function to read in files and arrange correctly ###
### updated to include all plotting functions #########
#######################################################

#' Graph Cerebral Blood Flow
#'
#' This function generates a perfusion map of cerebral blood flow.
#' It should be used for research and illustrative purposes only and should NOT be used to diagnose/treat any disease.
#' 
#' @param folder directory containing the set of CT images to be combined. Must be arranged in alphanumerical order according to the order in which they were collected
#' @param method character string indicating estimation technique to be used ("maxslope","sSVD","bSVD", or "oSVD")
#' @param aif character string indicating arterial input function to be used ("auto", "max", "gv", or "user")
#' @param t0 parameter for aif="gv". The time in which contrast enters. No default value. Required if aif="gv".
#' @param C0 parameter for aif="gv". Constant for gv function. Defaults to 1.
#' @param a parameter for aif="gv". Shape parameter for gv function. Defaults to 3.
#' @param b parameter for aif="gv". Rate parameter for gv function. Defaults to 1.5.
#' @param aif.user user-defined vector for arterial input function. Required if aif="user".
#' @param lag time interval in seconds for collection of CT data. Usually set to 1 or 2. Required.
#' 
#' @keywords ctp blood flow cbf perfusion
#' 
#' @examples
#' ctp_cbf(folder="C:/.../r001/loc12", method="maxslope", aif="gv", t0=10, lag=1) 
#' 
#' ctp_cbf(folder="C:/.../r083/loc3", method="oSVD", aif="gv", t0=10, C0=1.13, lag=1)
#' 
#' svec = c(rep(50,11),51,211,365,445,446,396,323,247,181,127,87,58,rep(50,4))
#' ctp_cbf(folder="C:/.../r083/loc3", method="oSVD", aif="user", aif.user=svec, lag=1)
#' 
#' @export
#' 

ctp_cbf = function(folder,method=c("maxslope","sSVD","bSVD","oSVD"), aif=c("auto", "max", "gv", "user"), t0, C0=1, a=3, b=1.5, aif.user, lag) {
  
  #setting color palette
  ctppalette = colorRampPalette(c("black", "purple4", "navy", "blue", "limegreen", "yellow2", "orange", "red"))
  
  filenames <- list.files(folder, pattern="*.dcm", full.names=TRUE)
  ldf <- lapply(filenames, readDICOMFile)
  lapply(ldf, names)
  
  ctarray = array(data=NA, dim=c(ncol(ldf[[1]]$img),nrow(ldf[[1]]$img),length(ldf)))
  for(i in 1:length(ldf)) {
    ctarray[, , i] = t(ldf[[i]]$img)
  }
  
  t = seq(from=0, to=(length(ldf)-1)*lag, by=lag)
  
  if(aif=="user") {
    aif.vec=aif.user
  }
  
  if(aif=="gv") {
    aif.f = function(t, t0, C0, a, b) {
      out = vector(mode="numeric", length=length(t))
      for(j in 1:length(t)) {
        i=t[j]
        if((i)<=t0)  val = 0
        if((i)>t0)   val = C0*(((i)-t0)^a)*exp((-((i)-t0))/b)
        k=val
        if(k<.5) k=.5 #baseine HU for tissue
        out[j] = k
      }
      return(out)
    }
    
    aif.vec = aif.f(t, t0, C0=C0, a=a, b=b)*100
  }
  
  if(aif=="max") {
    aif.vec = array(data=NA, dim=c(length(ldf)))
    for(i in 1:length(ldf)) {
      aif.vec[i] = max(ldf[[i]]$img)
    }
  }  
  
  if(aif=="auto") {
    #Arterial Input Block#
    deriv <- function(x, y) diff(y) / diff(x)
    middle_pts <- function(x) x[-1] - diff(x) / 2
    
    xdim = dim(ctarray)[1]
    ydim = dim(ctarray)[2]
    zdim = dim(ctarray)[3]
    
    caimage    = array(data=NA,dim=c(xdim,ydim))
    caimagen   = array(data=NA,dim=c(ceiling(xdim*ydim*.1),zdim))
    caimager   = array(data=NA,dim=c(ceiling(xdim*ydim*.1)))
    
    
    phantom2 = ctarray
    
    x=1:zdim
    id=order(x)
    
    for(i in 1:xdim) { 
      for(j in 1:ydim) {
        Ct = phantom2[i,j,]
        caimage[i,j] = sum(diff(x[id])*rollmean(Ct[id],2))
      }
    }
    
    locmat = which(caimage>=quantile(caimage,probs=0.9), arr.ind=T)
    
    caimage10  = array(data=NA,dim=c(nrow(locmat),zdim))
    
    g=1
    for(i in 1:nrow(locmat)) { 
      p = locmat[i,1]
      q = locmat[i,2]
      caimage10[g,] = phantom2[p,q,]
      g = g+1
    }
    
    
    for(i in 1:nrow(caimage10)) {
      xvec = 1:zdim
      id=order(xvec)
      deriv <- function(x, y) diff(y) / diff(x)
      middle_pts <- function(x) x[-1] - diff(x) / 2
      second_d <- deriv(middle_pts(xvec), deriv(xvec, caimage10[i,]))
      ss=second_d^2
      caimager[i] = trapz(ss,xvec) # this might not be the best, getting negative values
    }
    
    
    roughmat = which(caimager>=quantile(caimager,probs=0.25), arr.ind=T)
    
    caimager75 = array(data=NA,dim=c(length(roughmat),zdim))
    
    m=1
    for(i in 1:length(roughmat)) { 
      p = roughmat[i]
      caimager75[m,] = caimage10[p,]
      m = m+1
    }
    
    kclust = kmeans(caimager75[1:(nrow(caimager75)-1),1:length(t)], centers=5, nstart=10)
    
    clustnum = 0
    initmean = 100000000
    for(i in 1:5) {
      if(mean(kclust$centers[i,]) < initmean) {
        clustnum = i
        initmean = mean(kclust$centers[i,])
      }
    }
    
    kmat = cbind(caimager75[1:(nrow(caimager75)-1),], t(t(kclust$cluster)))
    
    kmatcn = which(kmat[,ncol(kmat)]==clustnum , arr.ind=T)
    
    kmat2 = array(data=NA,dim=c(length(kmatcn), length(t)))
    
    m=1
    for(i in 1:length(kmatcn)) { 
      p = kmatcn[i]
      kmat2[m,] = kmat[p,1:length(t)]
      m = m+1
    }
    
    
    kclust2 = kmeans(kmat2, centers=5, nstart=10)
    
    clustnum2 = 0
    initmean2 = 100000000
    for(i in 1:5) {
      if(mean(kclust2$centers[i,]) < initmean2) {
        clustnum2 = i
        initmean2 = mean(kclust2$centers[i,])
      }
    }
    
    kmat3 = cbind(kmat2, t(t(kclust2$cluster)))
    
    kmatcn2 = which(kmat3[,ncol(kmat3)]==clustnum2 , arr.ind=T)
    
    kmat4 = array(data=NA,dim=c(length(kmatcn2), length(t)))
    
    m=1
    for(i in 1:length(kmatcn2)) { 
      p = kmatcn2[i]
      kmat4[m,] = kmat3[p,1:length(t)]
      m = m+1
    }
    
    aif.vec = array(data=NA,dim=c(length(t)))
    
    for(i in 1:length(t)) {
      aif.vec[i] = mean(kmat4[,i])
    }
    
  }
  
  
  
  ######################
  
  #Max Slope block#
  if(method=="maxslope") {
    
    maxslope.cbf = function(Ca, inimage, lag=2) {
      xdim = dim(inimage)[1]
      ydim = dim(inimage)[2]
      zdim = dim(inimage)[3]
      
      outimage = array(data=NA,dim=c(xdim,ydim))
      
      for(i in 1:xdim) { 
        for(j in 1:ydim) {
          Ct = inimage[i,j,]
          outimage[i,j] = 60 * max(diff(Ct,lag=lag)) / max(Ca)
        }
      }
      
      #image.plot(outimage, col=ctppalette(100), axes=F, zlim=c(min(outimage),min(max(outimage),90)))
      brk = seq(0,75,5)
      image.plot(outimage, col=ctppalette(15), axes=F, zlim=c(0,75), breaks=brk, lab.breaks=names(brk), legend.lab="CBF")
      
    }
    
    maxslope.cbf(Ca=aif.vec, ctarray, lag=1)
  }
  
  #################
  
  #sSVD block#
  if(method=="sSVD") {     
    
    sSVD.cbf = function(Ca, inimage, lag=2) {
      xdim = dim(inimage)[1]
      ydim = dim(inimage)[2]
      zdim = dim(inimage)[3]
      
      camat = array(data=NA,dim=c(length(t),length(t)))
      for(i in 1:length(t)) {
        for(j in 1:length(t)) {
          if(j>i) camat[i,j] = 0
          if(i==j) camat[i,j] = Ca[1]
          if(i>j) camat[i,j] = Ca[(i-j)+1]
        }
      }
      
      cmat = camat
      cmat.svd = svd(cmat)
      u = (cmat.svd$u)
      v = (cmat.svd$v)
      d = (cmat.svd$d)
      psvd = .05*max(cmat.svd$d)
      
      sr = array(data=0,dim=c(length(t)))
      for(i in 1:length(sr)) {
        if(d[i]>psvd) sr[i] = 1/d[i]
      }
      s.inv = diag(sr)
      
      ca.inv = v %*% s.inv %*% t(u)
      
      outimage = array(data=NA,dim=c(xdim,ydim))
      
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          Ct = inimage[i,j,]
          rt.hat = (ca.inv %*% Ct) / lag
          cbf = max(rt.hat)
          outimage[i,j] = 22*cbf ## 10 gets it close
        }
      }
      
      #image.plot(outimage, col=ctppalette(100), axes=F, zlim=c(min(outimage),min(max(outimage),90)))
      brk = seq(0,75,5)
      image.plot(outimage, col=ctppalette(15), axes=F, zlim=c(0,75), breaks=brk, lab.breaks=names(brk), legend.lab="CBF")
    }
    
    sSVD.cbf(Ca=aif.vec, ctarray, lag=lag)
  }
  #################
  
  #bSVD block#
  if(method=="bSVD") {     
    
    bSVD.cbf = function(Ca, inimage, lag=2) {
      xdim = dim(inimage)[1]
      ydim = dim(inimage)[2]
      zdim = dim(inimage)[3]
      
      camat = array(data=NA,dim=c(length(t),length(t)))
      for(i in 1:length(t)) {
        for(j in 1:length(t)) {
          if(j>i) camat[i,j] = 0
          if(i==j) camat[i,j] = Ca[1]
          if(i>j) camat[i,j] = Ca[(i-j)+1]
        }
      }
      
      B = array(data=0,dim=c(length(t),length(t)))
      for(i in 1:length(t)) {
        for(j in 1:length(t)) {
          if(j>i) B[i,j] = Ca[length(t)-(j-i)+1]
        }
      }
      
      cmat = rbind(cbind(camat,B),cbind(B,camat))
      cmat.svd = svd(cmat)
      u = (cmat.svd$u)
      v = (cmat.svd$v)
      d = (cmat.svd$d)
      psvd = .05*max(cmat.svd$d)
      
      sr = array(data=0,dim=c(length(cmat[1,])))
      for(i in 1:length(sr)) {
        if(d[i]>psvd) sr[i] = 1/d[i]
      }
      s.inv = diag(sr)
      
      ca.inv = v %*% s.inv %*% t(u)
      
      outimage = array(data=NA,dim=c(xdim,ydim))
      
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          Ct = inimage[i,j,]
          Ctb = rep(0,length(cmat[1,]))
          for(k in 1:length(Ct)) {Ctb[k]=Ct[k]} ## padding zeros at the end
          rt.hat = (ca.inv %*% Ctb) / lag
          cbf = max(rt.hat)
          outimage[i,j] = 80*cbf ## 40 gets it close
        }
      }
      
      #image.plot(outimage, col=ctppalette(100), axes=F, zlim=c(min(outimage),min(max(outimage),90)))
      brk = seq(0,75,5)
      image.plot(outimage, col=ctppalette(15), axes=F, zlim=c(0,75), breaks=brk, lab.breaks=names(brk), legend.lab="CBF")
    }
    
    bSVD.cbf(Ca=aif.vec, ctarray, lag=lag)
  }
  #################
  
  #oSVD block#
  if(method=="oSVD") {    
    
    oSVD.cbf = function(Ca, inimage, lag=2, oi) {
      xdim = dim(inimage)[1]
      ydim = dim(inimage)[2]
      zdim = dim(inimage)[3]
      
      camat = array(data=NA,dim=c(length(t),length(t)))
      for(i in 1:length(t)) {
        for(j in 1:length(t)) {
          if(j>i) camat[i,j] = 0
          if(i==j) camat[i,j] = Ca[1]
          if(i>j) camat[i,j] = Ca[(i-j)+1]
        }
      }
      
      B = array(data=0,dim=c(length(t),length(t)))
      for(i in 1:length(t)) {
        for(j in 1:length(t)) {
          if(j>i) B[i,j] = Ca[length(t)-(j-i)+1]
        }
      }
      
      cmat = rbind(cbind(camat,B),cbind(B,camat))
      cmat.svd = svd(cmat)
      u = (cmat.svd$u)
      v = (cmat.svd$v)
      d = (cmat.svd$d)
      psvd = .05*max(cmat.svd$d)
      
      sr = array(data=0,dim=c(length(cmat[1,])))
      for(i in 1:length(sr)) {
        if(d[i]>psvd) sr[i] = 1/d[i]
      }
      s.inv = diag(sr)
      
      ca.inv = v %*% s.inv %*% t(u)
      
      outimage = array(data=NA,dim=c(xdim,ydim))
      
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          Ct = inimage[i,j,]
          #Ct[Ct>2000] = NA
          Ctb = rep(0,length(cmat[1,]))
          for(k in 1:length(Ct)) {Ctb[k]=Ct[k]} ## padding zeros at the end
          rt.hat = (ca.inv %*% Ctb) / lag
          sumabsf = 0
          for(v in 3:(length(cmat[1,])-1)) {sumabsf = sumabsf + abs(rt.hat[v] - 2*rt.hat[v-1] + rt.hat[v-2])}
          cbf = (1/length(cmat[1,])) * (1/oi) * (sumabsf)
          outimage[i,j] = 60*cbf 
        }
      }
      #image.plot(outimage, col=ctppalette(100), axes=F, zlim=c(min(outimage),min(max(outimage),90)))
      brk = seq(0,75,5)
      image.plot(outimage, col=ctppalette(15), axes=F, zlim=c(0,75), breaks=brk, lab.breaks=names(brk), legend.lab="CBF")
    }
    
    oSVD.cbf(Ca=aif.vec, ctarray, oi=0.2, lag=lag)
  }
  #################
  
  
}
#######################################################
#######################################################
