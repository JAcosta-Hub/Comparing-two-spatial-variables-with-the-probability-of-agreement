##################################################################
##--------------------------------------------------------------##
## Spatiotemporal Probability of Agreement for simulations 

#h: spatial lag
#u: temporal lag
#C: value of the constant c (CAD)
#trend_t: type of temporal trend (constant and linear permitted)
#nugget: value of the nugget parameter
Psi_est_st <- function(h=0, u=0, C=1, param.est, corrmodel, trend_t="linear", nugget=0, ...){
  Psi_c=c()
  for(j in 1:nrow(param.est)){
    param=c(param.est[j, ], nugget=nugget)
    Psi_c[j]=Psi_ST(h = h, u = u, C = C,
                    paramcorr = param, corrmodel = corrmodel, trend_t = trend_t)
  }
  return(Psi_c)
}

##--------------------------------------------------------------##
## Spatiotemporal Probability of Agreement 

#h: spatial lag
#u: temporal lag
#C: value of the constant c (CAD)
#paramcorr: values of the covariance model parameters
#corrmodel: covariance models (exponential separable and Iacocesare permitted)
#trend_t: type of temporal trend (constant and linear permitted)
Psi_ST <- function(h, u, C=2, paramcorr, corrmodel, trend_t="linear",...){
  if(trend_t=="const"){muD=0}
  if(trend_t=="linear"){muD=u*paramcorr$mean1}
  sigmaDh=sqrt(sigmaD_ST(h, u, paramcorr, corrmodel=corrmodel))
  x1=(C-muD)/sigmaDh
  x2=-(C+muD)/sigmaDh
  if(h==0&u==0&C==0){res=0}else{res=pnorm(x1)-pnorm(x2)}
  return(res)
}


##-------------------------------------------------------##
## Spatiotemporal covariance functions of differences

#h: spatial lag
#u: temporal lag
#paramcorr: values of the covariance model parameters
#corrmodel: covariance models (exponential separable and Iacocesare permitted)
sigmaD_ST <- function(h, u, paramcorr, corrmodel){
  C_00=paramcorr$sill+paramcorr$nugget
  hs=h/paramcorr$scale_s
  ut=u/paramcorr$scale_t
  if(corrmodel=="Exp_Exp"){ C_hu=paramcorr$sill*exp(-hs)*exp(-ut) }
  if(corrmodel=="Iacocesare"){C_hu=paramcorr$sill*(1+hs^paramcorr$power_s+ut^paramcorr$power_t)^(-paramcorr$power2)   }
  sigmaD_ST=2*(C_00-C_hu)
  return(sigmaD_ST)
}


#######################################################
##---------------------------------------------------##
## Gcc funstions

#ima: list of all images in RGB format to be processed, i.e. each image is an array of dimension c(n,m,3)
#trend_t: type of temporal trend (constant and linear permitted)
gcc_st <- function(ima, trend_t="linear"){
  NT=length(ima)
  aux1=dim(ima[[1]])
  nx=aux1[1]
  ny=aux1[2]
  
  times <- seq(1, NT, 1)
  coords=expand.grid(x=0:(nx-1), y=0:(ny-1))
  if(trend_t=="linear"){ X=cbind(rep(1,nx*ny*NT), rep(times, each=nx*ny)) }else{ X=NULL }
  
  aux2=matrix(NA, nr=NT, nc=nx*ny)
  for(j in 1:NT){
    aux3=ima[[j]][,,2]/(ima[[j]][,,1]+ima[[j]][,,2]+ima[[j]][,,3])
    aux2[j, ] <- c(aux3)
  }
  return(list(data=aux2, coords=coords, times=times, X=X))
}

##-------------------------------------##
## Rasterization

#x: list of all images in RGB format to be rasterized, i.e. each image is an array of dimension c(n,m,3)
#px: number of pixels to be grouped in the form px^2
#raster_name: name with which the raster images will be saved in the PC, if it is NULL it does not save these images.
rasterization <- function(x, px=10, raster_name=NULL){
  if(!is.list(x)){x=list(x)}
  aux=trunc(dim(x[[1]])/px)
  n1=aux[1]
  n2=aux[2]
  raster_image <- list()
  for(k in 1:length(x)){
    aux=array(NA, dim=c(n1,n2,3))
    for(i in 1:n1){
      for(j in 1:n2){
        aux[i,j,1]<- mean(x[[k]][(1:px)+(i-1)*px, 1:px+(j-1)*px, 1])
        aux[i,j,2]<- mean(x[[k]][(1:px)+(i-1)*px, 1:px+(j-1)*px, 2])
        aux[i,j,3]<- mean(x[[k]][(1:px)+(i-1)*px, 1:px+(j-1)*px, 3])
      }
    }
    raster_image[[k]] <- aux
  }
  if(!is.null(raster_name)){
    for(i in 1:length(raster_name)){
      jpeg(filename=raster_name[i],
           width = n2, height = n1, quality = 100)
      OpenImageR::imageShow(raster_image[[i]])
      dev.off()
    }
  }
  return(raster_image)
}


################################################################
##----------------------------------------------##
## Spatiotemporal model parameter estimation

#sim: is a list containing: the data, the coordinates, the model, the parameters with which it was simulated. 
#name.fixed: name of the parameters that will be fixed, i.e., that will not be estimated
#name.start: name of the parameters that will be estimated.
est_param_st<-function(sim, name.fixed, name.start){
  NS=ncol(sim$data.sim[[1]][[1]])
  NT=nrow(sim$data.sim[[1]][[1]])
  
  # Fixed parameters and starting value for the estimated parameters
  fixed<-list()
  start<-list()
  for(k in 1:length(sim$param.sim)){    
    aux=data.frame(t(unlist(sim$param.sim[[k]])))
    fixed[[k]] <- param.form2(aux[name.fixed])
    start[[k]] <- param.form2(aux[name.start])
  }
  
  # Maximum composite-likelihood fitting of the RF:
  fit_ST <- list()
  for(k in 1:length(start) ){
    fit.aux <- list()
    for(j in 1: length(sim[[1]][[1]]) ){
      t <- proc.time()
      print(paste("caso", k, "iteración", j))
      fit.aux[[j]] <- GeoFit(data=sim$data.sim[[k]][[j]],
                             coordx=sim$coords.sim, 
                             coordt=sim$times.sim,
                             corrmodel=sim$corrmodel, 
                             maxtime=2, maxdist=start[[k]]$scale_s/10,
                             likelihood="Marginal", type="Pairwise",
                             start=start[[k]], fixed=fixed[[k]], 
                             X=sim$X)
      t1 <- proc.time()-t
      print(t1)
    }
    fit_ST[[k]] <- fit.aux
  }
  return(fit_ST)
}

################################################################
##------------------------------------------------------------##
## Plot of simulations of the spatiotemporal process

# data.sim: list of simulations 
# name.sim: name under which the graphic will be saved 
# NS: spatial dimension, of the form NS^2
# opcion: type of image to be saved, option=1, is the normal image, option=2, corresponds to a smoothed version by kriging
plot_sim <- function(data.sim, name.sim, NS=50, opcion=1){
  m1=min(data.sim[[1]][[1]])
  m2=max(data.sim[[1]][[1]])
  if(opcion==1){
    png(name.sim, width=1024, height=768, res=130)
    par(mfrow=c(2,3), mar=c(2,2,2,1))
    image(matrix(data.sim[[1]][[1]][1,],NS,NS), xaxt="n", yaxt="n", main=expression(t==1), zlim=c(m1,m2), col=terrain.colors(64))
    image(matrix(data.sim[[1]][[1]][2,],NS,NS), xaxt="n", yaxt="n", main=expression(t==2), zlim=c(m1,m2), col=terrain.colors(64))
    image(matrix(data.sim[[1]][[1]][3,],NS,NS), xaxt="n", yaxt="n", main=expression(t==3), zlim=c(m1,m2), col=terrain.colors(64))
    image(matrix(data.sim[[1]][[1]][4,],NS,NS), xaxt="n", yaxt="n", main=expression(t==4), zlim=c(m1,m2), col=terrain.colors(64))
    image(matrix(data.sim[[1]][[1]][5,],NS,NS), xaxt="n", yaxt="n", main=expression(t==5), zlim=c(m1,m2), col=terrain.colors(64))
    image(matrix(data.sim[[1]][[1]][6,],NS,NS), xaxt="n", yaxt="n", main=expression(t==6), zlim=c(m1,m2), col=terrain.colors(64))
    dev.off()
  }
  if(opcion==2){
    ks=krige.conv(s100, locations=expand.grid(1:NS, 1:NS), krige=krige.control(cov.pars=c(1, .25)) )
    png(name.sim, width=1024, height=768, res=130)
    par(mfrow=c(2,3), mar=c(2,2,2,1))
    contour(ks, values=data.sim[[1]][[1]][1,], main=expression(t==1), color=terrain.colors, filled=TRUE, zlim=c(m1,m2))
    contour(ks, values=data.sim[[1]][[1]][2,], main=expression(t==2), color=terrain.colors, filled=TRUE, zlim=c(m1,m2))
    contour(ks, values=data.sim[[1]][[1]][3,], main=expression(t==3), color=terrain.colors, filled=TRUE, zlim=c(m1,m2))
    contour(ks, values=data.sim[[1]][[1]][4,], main=expression(t==4), color=terrain.colors, filled=TRUE, zlim=c(m1,m2))
    contour(ks, values=data.sim[[1]][[1]][5,], main=expression(t==5), color=terrain.colors, filled=TRUE, zlim=c(m1,m2))
    contour(ks, values=data.sim[[1]][[1]][6,], main=expression(t==6), color=terrain.colors, filled=TRUE, zlim=c(m1,m2))
    dev.off()
  }
}

##-------------------------------------------##
## simulation of space-time random field

#NS: Size of the square grid of the form NS^2
#NT: Size of the time dimension
#M: Number of simulations
#paramcorr: values of the covariance model parameters
#corrmodel: covariance models (exponential separable and Iacocesare permitted)
#trend_t: type of temporal trend (constant and linear permitted)
#name.sim: Name under which the simulations will be saved on the disk
#plot.it: If true instead of saving the simulations it saves a realization. NT must >=6 and it is suggested to use with M=1.
sim_st <- function(NS, NT, M, param.sim, corrmodel, trend_t, name.sim, plot.it=FALSE){
 ## Spatial coordinates - Temporal coordinates -  Time trend
 coords.sim=expand.grid(x=0:(NS-1), y=0:(NS-1))
 times.sim <- seq(1, NT, 1)
 if(trend_t=="linear"){ X=cbind(rep(1,NS^2*NT), rep(times.sim, each=NS^2)) }else{ X=NULL }

 ## Simulations
 data.sim=list()
 for(k in 1:length(param.sim)){
  data_aux=list()
  for(j in 1:M){
    t <- proc.time()
    print(paste("caso", k, "iteración", j))
    data_aux[[j]] <- GeoSim(coordx=coords.sim, coordt=times.sim, corrmodel=corrmodel, 
                            param=param.sim[[k]], X=X)$data
    t1 <- proc.time()-t
    print(t1)
    }
  data.sim[[k]] <- data_aux 
 }
 
 if(plot.it){
   ## Save maps of simulations
   png(name.sim, width=1024, height=960, res=130)
   par(mfrow=c(2,3), mar=c(2,2,2,1))
   m1=min(data.sim[[1]][[1]])
   m2=max(data.sim[[1]][[1]])
   image(matrix(data.sim[[1]][[1]][1,],NS,NS), xaxt="n", yaxt="n", main=expression(t==1), zlim=c(m1,m2), col=terrain.colors(64))
   image(matrix(data.sim[[1]][[1]][2,],NS,NS), xaxt="n", yaxt="n", main=expression(t==2), zlim=c(m1,m2), col=terrain.colors(64))
   image(matrix(data.sim[[1]][[1]][3,],NS,NS), xaxt="n", yaxt="n", main=expression(t==3), zlim=c(m1,m2), col=terrain.colors(64))
   image(matrix(data.sim[[1]][[1]][4,],NS,NS), xaxt="n", yaxt="n", main=expression(t==4), zlim=c(m1,m2), col=terrain.colors(64))
   image(matrix(data.sim[[1]][[1]][5,],NS,NS), xaxt="n", yaxt="n", main=expression(t==5), zlim=c(m1,m2), col=terrain.colors(64))
   image(matrix(data.sim[[1]][[1]][6,],NS,NS), xaxt="n", yaxt="n", main=expression(t==6), zlim=c(m1,m2), col=terrain.colors(64))
   dev.off()
   return(data.sim)
 }else{
  ## Save simulations
  save(data.sim, coords.sim, times.sim, param.sim, corrmodel, X, file = name.sim) 
  return(list(data.sim=data.sim, 
             coords.sim=coords.sim, 
             times.sim=times.sim, 
             param.sim=param.sim,
             corrmodel=corrmodel, 
             X=X) ) 
 } 
}

##-------------------------------------------------------##
## Format Initial parameters
#param: Parameters in data.frame formt
param.form <- function(param){
  name<-names(param)
  param.name=paste0("param",1:nrow(param))
  param.sim=list()
  print(nrow(param))
  for(i in 1:nrow(param)){
    aux1=list()
    for(j in 1:ncol(param)){
      aux1[[name[j]]] <- param[i, j]
    }
    param.sim[[param.name[i]]] <- aux1
  }
  print(param.sim)  
  return(param.sim)
}

param.form2 <- function(param){
  name<-names(param)
  param=as.matrix(param)
  print(dim(param))
  print(name)
  param.sim=list()
    for(j in 1:ncol(param)){
      param.sim[[name[j]]] <- as.numeric(param[1, j])
    }
  return(param.sim)
}


##--------------------------------------##
## View Estimated Parameters
#x: Parameters estimated using est_param_st function.
ver_param <- function(x){
  m=length(x)
  param.est=c()
  for(j in 1:m){
    param.est=rbind(param.est, unlist(x[[j]]$param) )
  }
  param.est <- data.frame(param.est)
  #   colnames(param.est) <-  names(x[[1]]$param)
  return(param.est)
}

ver_param_filtro <- function(x, cond, param_fil){
  x_fil <- list()
  for(k in 1:length(x)){ 
     if(!is.data.frame(x[[k]])){
       y=ver_param(x[[k]])[,param_fil ]
       fil=filtro(x=y, cond=cond)
       x_fil[[k]] <- ver_param(x[[k]])[fil==1,]
     }else{
       y=x[[k]][,param_fil]
       fil=filtro(x=y, cond=cond)
       x_fil[[k]] <- x[[k]][fil==1,]
     }
  }
  return(x_fil)
}

##----------------------------------------##
## Filter out unusual estimators  
filtro <- function(x, cond)
{
  m=length(cond)
  st=substr(cond,1,1)
  num=as.numeric(substring(cond,2))
  
  fil=matrix(NA, ncol=m, nrow=length(x))
  for(j in 1:m){
    if(st[j]=="<"){ fil[,j] <- !(x<num[j]) }
    if(st[j]==">"){ fil[,j] <- !(x>num[j]) }
  } 
  return(apply(fil,1,prod))
}



###############################################
##-------------------------------------------##
## Probabilidad de Agreement Bivariado

Psi_est_biv <- function(h=0, C=1, est, ...){
  Psi_hc=c()
  for(j in 1:nrow(est)){
    Psi_hc[j]=Psi(h=h, C=C,
                  muD=(est$mean_1-est$mean_2)[j], 
                  sigma2=c(est$sill_1[j], est$sill_2[j]), 
                  phi=c(est$scale_1[j], est$scale_2[j], est$scale_12[j]), 
                  rhoxy=est$pcol[j])
  }
  return(Psi_hc)
}

##################################
## Probability of Agreement
Psi <- function(h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.2), nu=c(0.5,0.5,0.5), rhoxy=0, ...){
  sigmaDh=sqrt(sigmaD(h, sigma2=sigma2, phi=phi, nu=nu, rhoxy=rhoxy))
  x1=(C-muD)/sigmaDh
  x2=-(C+muD)/sigmaDh
  return(pnorm(x1)-pnorm(x2))
}


#########################################
## Covariance of bivariate Difference
sigmaD <- function(h, sigma2, phi, nu, rhoxy){
  Cxy=rhoxy*sqrt(sigma2[1]*sigma2[2])*cov.spatial(h, cov.model="matern", cov.pars=c(1, phi[3]), kappa=nu[3])
  sigmaD=sum(sigma2[1:2])-2*Cxy
  return(sigmaD)
}

#################################################
## Bivariate random number simulation
sim_bivariado <- function(NN, M, param, name.sim){
  ##--------------------------------------##  
  ## Estructura del Modelo    
  coords.sim=expand.grid(coordx=0:(NN-1) , coordy=0:(NN-1) )
  model="Gaussian" 
  corrmodel="Bi_Matern"
  
  ##--------------------------------------##  
  ## Simulaciones 
  data.sim=list()
  for( k in 1:length(param)){
    ss.aux=list()
    t1=list()
    for(j in 1:M){
      t <- proc.time()
      print(paste("caso", k, "iteracion", j))
      ss.aux[[j]] = GeoSim(coordx=coords.sim, corrmodel=corrmodel, model=model, 
                           param=param[[k]])$data
      t1[[j]] <- proc.time()-t
      print(t1[[j]]) }
    data.sim[[k]] <- ss.aux }
  
  ##--------------------------------------##
  ## Guardar Simulaciones
  save(data.sim, coords.sim, param.sim, file = name.sim) 
  
  return(list(data.sim=data.sim, 
              coords.sim=coords.sim, 
              param.sim=param.sim))
}     

