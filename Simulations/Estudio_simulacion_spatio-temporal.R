##########################################
##--------------------------------------##
## Librerias y funciones

library(GeoModels)
library(geoR)
library(xtable)
source("Functions_prob_agreement.R")

 ##########################################
 ##--------------------------------------##
 ## Initial parameters

 ## Separable Model
 CorrParam("Exp_Exp")
 param0.exp_exp=expand.grid(nugget=0, mean=0.5, mean1=c(-0.1, 0.1),
                            sill=0.1, scale_s=20/log(1/0.05), scale_t=1 )

 param0.exp_exp
 param.sim.exp_exp=param.form(param0.exp_exp)
 length(param.sim.exp_exp)
 
 ## Nonseparable Model
 CorrParam("Iacocesare")
 param0.iacocesare=expand.grid(nugget=0, mean=0.5, mean1=c(-0.1, 0.1),
                            sill=0.1, scale_s=20/log(1/0.05), scale_t=1,
                            power2=2, power_s=1, power_t=1)
 param0.iacocesare
 param.sim.iacocesare=param.form(param0.iacocesare)
 length(param.sim.iacocesare)
 
 ##########################################
 ##--------------------------------------##
 ## Graphic of a simulation
 
 sim_st_exp_exp=sim_st(NS=50, NT=6, M=1, param.sim=list(param.sim.exp_exp[[2]]),
        corrmodel="Exp_Exp", trend_t="linear", plot.it=TRUE,
        name.sim="relisation_ST_NS50_exp_exp.png")

 sim_st_iacocesare=sim_st(NS=50, NT=6, M=1, param.sim=list(param.sim.iacocesare[[2]]),
        corrmodel="Iacocesare", trend_t="linear", plot.it=TRUE,
        name.sim="relisation_ST_NS50_iacocesare.png")
 
 # If you already have the simulation you can plot directly
 plot_sim(sim_st_exp_exp, "relisation_ST_NS50_exp_exp.png")
 plot_sim(sim_st_iacocesare, "relisation_ST_NS50_iacocesare.png")
 

 ##########################################
 ##--------------------------------------##
 ## Simulations
 #If you run this part it can be slow
 sim_NS20_T10_separable=sim_st(NS=20, NT=10, M=250, param.sim=param.sim.exp_exp,
                         corrmodel="Exp_Exp", trend_t="linear", 
                         name.sim="Sim_NS20_T10_separable.RData")
 
 sim_NS20_T10_nonsepar=sim_st(NS=20, NT=10, M=250, param.sim=param.sim.iacocesare,
                            corrmodel="Iacocesare", trend_t="linear", 
                            name.sim="Sim_NS20_T10_nonsepar.RData")
 

 sim_NS50_T10_separable=sim_st(NS=50, NT=10, M=250, param.sim=param.sim.exp_exp,
                               corrmodel="Exp_Exp", trend_t="linear", 
                               name.sim="Sim_NS50_T10_separable.RData")
 
 sim_NS50_T10_nonsepar=sim_st(NS=50, NT=10, M=250, param.sim=param.sim.iacocesare,
                              corrmodel="Iacocesare", trend_t="linear", 
                              name.sim="Sim_NS50_T10_nonsepar.RData")
 
 ##############################################
 ##------------------------------------------##
 ## Parameter Estimation: Separable Case
 # If you do not execute the above, it is recommended to load the fields already simulated to verify results.
 # load("Sim_NS20_T10_separable.RData")
 # load("Sim_NS50_T10_separable.RData")
 fit_NS20_T10_separable <- est_param_st(sim=sim_NS20_T10_separable, 
                                       name.fixed=names(param0.exp_exp)[1],
                                       name.start=names(param0.exp_exp)[-1])
 
 fit_NS50_T10_separable <- est_param_st(sim=sim_NS50_T10_separable, 
                                        name.fixed=names(param0.exp_exp)[1],
                                        name.start=names(param0.exp_exp)[-1])
 
 ## Visualization of parameters estimated
 par(mfrow=c(2,2))
 boxplot(ver_param(fit_NS20_T10_separable[[1]]), main="NS=20 & a1<0" )
 points( unlist(param0.exp_exp[2,names(fit_NS20_T10_separable[[2]][[1]]$param)]), pch=20, col=2, cex=2)
 
 boxplot(ver_param(fit_NS20_T10_separable[[2]]), main="NS=20 & a1>0" )
 points( unlist(param0.exp_exp[2,names(fit_NS20_T10_separable[[2]][[1]]$param)]), pch=20, col=2, cex=2)
 
 boxplot(ver_param(fit_NS50_T10_separable[[1]]), main="NS=50 & a1<0" )
 points( unlist(param0.exp_exp[1,names(fit_NS50_T10_separable[[1]][[1]]$param)]), pch=20, col=2, cex=2)
 
 boxplot(ver_param(fit_NS50_T10_separable[[2]]), main="NS=50 & a1>0" )
 points( unlist(param0.exp_exp[2,names(fit_NS50_T10_separable[[2]][[1]]$param)]), pch=20, col=2, cex=2)
 
 ##----------------------------------------##
 ## Eliminate rare cases
 
 ## condition to filter out rare cases in phi_s
 NN=20
 cond1=c("<0", paste0(">", round(NN*sqrt(2)/log(1/0.05)*2,1)) )
 
 fit2_NS20_T10_separable <- ver_param_filtro(x=fit_NS20_T10_separable, cond=cond1, param_fil="scale_s")
 fit2_NS50_T10_separable <- ver_param_filtro(x=fit_NS50_T10_separable, cond=cond1, param_fil="scale_s")
 
 ## Visualization
 par(mfrow=c(2,2))
 boxplot(fit2_NS20_T10_separable[[1]], main="NS=20 & a1<0" )
 points( unlist(param0.exp_exp[1,names(fit2_NS20_T10_separable[[1]])]), pch=20, col=2, cex=2)
 
 boxplot(fit2_NS20_T10_separable[[2]], main="NS=20 & a1>0"  )
 points( unlist(param0.exp_exp[2,names(fit2_NS20_T10_separable[[2]])]), pch=20, col=2, cex=2)
 
 boxplot(fit2_NS50_T10_separable[[1]], main="NS=50 & a1<0"  )
 points( unlist(param0.exp_exp[1,names(fit2_NS50_T10_separable[[1]])]), pch=20, col=2, cex=2)
 
 boxplot(fit2_NS50_T10_separable[[2]], main="NS=50 & a1<0"  )
 points( unlist(param0.exp_exp[2,names(fit2_NS50_T10_separable[[2]])]), pch=20, col=2, cex=2)
 
 
 #################################################
 ##---------------------------------------------##
 ## Estimación de Parámetros: Caso No Separable
 # If you do not execute the above, it is recommended to load the fields already simulated to verify results.
 # load("Sim_NS20_T10_nonsepar.RData")
 # load("Sim_NS50_T10_nonsepar.RData")
 fit_NS20_T10_nonsepar <- est_param_st(sim=sim_NS20_T10_nonsepar, 
                                       name.fixed=names(param0.iacocesare)[1],
                                       name.start=names(param0.iacocesare)[-1])
 
 fit_NS50_T10_nonsepar <- est_param_st(sim=sim_NS50_T10_nonsepar, 
                                       name.fixed=names(param0.iacocesare)[1],
                                       name.start=names(param0.iacocesare)[-1])
 
 ## Visualization of parameters estimated
 par(mfrow=c(2,2))
 boxplot(ver_param(fit_NS20_T10_nonsepar[[1]]), main="NS=20 & a1<0" )
 points( unlist(param0.iacocesare[1,names(fit_NS20_T10_nonsepar[[1]][[1]]$param)]), pch=20, col=2, cex=2)
 
 boxplot(ver_param(fit_NS20_T10_nonsepar[[2]]), main="NS=20 & a1>0"  )
 points( unlist(param0.iacocesare[2,names(fit_NS50_T10_nonsepar[[2]][[1]]$param)]), pch=20, col=2, cex=2)
 
 boxplot(ver_param(fit_NS50_T10_nonsepar[[1]]), main="NS=50 & a1<0" )
 points( unlist(param0.iacocesare[1,names(fit_NS50_T10_nonsepar[[1]][[1]]$param)]), pch=20, col=2, cex=2)
 
 boxplot(ver_param(fit_NS50_T10_nonsepar[[2]]), main="NS=50 & a1>0"  )
 points( unlist(param0.iacocesare[2,names(fit_NS50_T10_nonsepar[[2]][[1]]$param)]), pch=20, col=2, cex=2)
 
 ## Eliminate rare cases
 # The same condition applies as in the separable case

 fit2_NS20_T10_nonsepar <- ver_param_filtro(x=fit_NS20_T10_nonsepar, cond=cond1, param_fil="scale_s") 
 fit2_NS50_T10_nonsepar <- ver_param_filtro(x=fit_NS50_T10_nonsepar, cond=cond1, param_fil="scale_s")
 
 ## Visualization of parameters estimated
 par(mfrow=c(2,2))
 boxplot(fit2_NS20_T10_nonsepar[[1]], main="NS=20 & a1<0"  )
 points( unlist(param0.iacocesare[1, names(fit2_NS50_T10_nonsepar[[1]])]), pch=20, col=2, cex=2)
 
 boxplot(fit2_NS20_T10_nonsepar[[2]], main="NS=20 & a1>0"  )
 points( unlist(param0.iacocesare[2,names(fit2_NS50_T10_nonsepar[[2]])]), pch=20, col=2, cex=2)
 
 boxplot(fit2_NS50_T10_nonsepar[[1]], main="NS=50 & a1<0"  )
 points( unlist(param0.iacocesare[1, names(fit2_NS50_T10_nonsepar[[1]])]), pch=20, col=2, cex=2)
 
 boxplot(fit2_NS50_T10_nonsepar[[2]], main="NS=50 & a1>0"  )
 points( unlist(param0.iacocesare[2,names(fit2_NS50_T10_nonsepar[[2]])]), pch=20, col=2, cex=2)
 

 ## condition to filter out rare cases in power_s and power_t
 
 cond2=c("<0", paste0(">", 10) )
 fit3_NS20_T10_nonsepar <- ver_param_filtro(x=fit2_NS20_T10_nonsepar, cond=cond2, param_fil="power_s") 
 fit3_NS50_T10_nonsepar <- ver_param_filtro(x=fit2_NS50_T10_nonsepar, cond=cond2, param_fil="power_s") 
 
 apply(fit3_NS20_T10_nonsepar[[1]], 2, min)
 apply(fit3_NS20_T10_nonsepar[[1]], 2, max)
 
 fit4_NS20_T10_nonsepar <- ver_param_filtro(x=fit3_NS20_T10_nonsepar, cond=cond2, param_fil="power_t") 
 fit4_NS50_T10_nonsepar <- ver_param_filtro(x=fit3_NS50_T10_nonsepar, cond=cond2, param_fil="power_t") 
 
 apply(fit4_NS20_T10_nonsepar[[1]], 2, min)
 apply(fit4_NS20_T10_nonsepar[[1]], 2, max)
 
 nrow(fit4_NS20_T10_nonsepar[[1]])/250*100
 nrow(fit4_NS20_T10_nonsepar[[2]])/250*100

 nrow(fit4_NS50_T10_nonsepar[[1]])/250*100
 nrow(fit4_NS50_T10_nonsepar[[2]])/250*100
 
  
 ## Visualization of parameters estimated
 par(mfrow=c(2,2))
 boxplot(fit4_NS20_T10_nonsepar[[1]], main="NS=20 & a1<0"  )
 points( unlist(param0.iacocesare[1, names(fit2_NS50_T10_nonsepar[[1]])]), pch=20, col=2, cex=2)
 
 boxplot(fit4_NS20_T10_nonsepar[[2]], main="NS=20 & a1>0"  )
 points( unlist(param0.iacocesare[2,names(fit2_NS50_T10_nonsepar[[2]])]), pch=20, col=2, cex=2)
 
 boxplot(fit4_NS50_T10_nonsepar[[1]], main="NS=50 & a1<0"  )
 points( unlist(param0.iacocesare[1, names(fit2_NS50_T10_nonsepar[[1]])]), pch=20, col=2, cex=2)
 
 boxplot(fit4_NS50_T10_nonsepar[[2]], main="NS=50 & a1>0"  )
 points( unlist(param0.iacocesare[2,names(fit2_NS50_T10_nonsepar[[2]])]), pch=20, col=2, cex=2)
 
 
 #######################################################
 ##---------------------------------------------------##
 ## Results of Estimation: Separable case

 # Percent valid negative slope
 nrow(fit2_NS20_T10_separable[[1]])/250*100
 nrow(fit2_NS50_T10_separable[[1]])/250*100

 # Percent valid positive slope 
 nrow(fit2_NS20_T10_separable[[2]])/250*100
 nrow(fit2_NS50_T10_separable[[2]])/250*100
 
 # Estimates
 print(booktab=TRUE, xtable(rbind(
    true1=unlist(param0.exp_exp[1,names(fit2_NS20_T10_separable[[1]])]),
    mean1.20=apply(fit2_NS20_T10_separable[[1]], 2, mean),
    sd1.20=apply(fit2_NS20_T10_separable[[1]], 2, sd),
    mean1.50=apply(fit2_NS50_T10_separable[[1]], 2, mean),
    sd1.50=apply(fit2_NS50_T10_separable[[1]], 2, sd),
    true2=unlist(param0.exp_exp[2,names(fit2_NS20_T10_separable[[1]])]),
    mean2.20=apply(fit2_NS20_T10_separable[[2]], 2, mean),
    sd2.20=apply(fit2_NS20_T10_separable[[2]], 2, sd),
    mean2.50=apply(fit2_NS50_T10_separable[[2]], 2, mean),
    sd2.50=apply(fit2_NS50_T10_separable[[2]], 2, sd)), digits=3 ))
 
 ##---------------------------------------------------##
 ## Results of Estimation: Nonseparable case

 # Percent valid negative slope
 nrow(fit4_NS20_T10_nonsepar[[1]])/250*100
 nrow(fit4_NS50_T10_nonsepar[[1]])/250*100
 
 # Percent valid positive slope
 nrow(fit4_NS20_T10_nonsepar[[2]])/250*100
 nrow(fit4_NS50_T10_nonsepar[[2]])/250*100
 
 # Estimates
 print(booktab=TRUE, xtable(rbind(
    true1=unlist(param0.iacocesare[1,names(fit2_NS50_T10_nonsepar[[1]])]),
    mean1.20=apply(fit4_NS20_T10_nonsepar[[1]], 2, mean),
    sd1.20=apply(fit4_NS20_T10_nonsepar[[1]], 2, sd),
    mean1.50=apply(fit4_NS50_T10_nonsepar[[1]], 2, mean),
    sd1.50=apply(fit4_NS50_T10_nonsepar[[1]], 2, sd),
    true2=unlist(param0.iacocesare[2,names(fit2_NS50_T10_nonsepar[[1]])]),
    mean2.20=apply(fit4_NS20_T10_nonsepar[[2]], 2, mean),
    sd2.20=apply(fit4_NS20_T10_nonsepar[[2]], 2, sd),
    mean2.50=apply(fit4_NS50_T10_nonsepar[[2]], 2, mean),
    sd2.50=apply(fit4_NS50_T10_nonsepar[[2]], 2, sd))[,names(fit2_NS50_T10_nonsepar[[1]])[c(1,2,6:8,3:5)]], digits=3 ))
 
 
 #######################################################
 ##---------------------------------------------------##
 ## Probability of Agreement: Variation of h, u and c
 
 h2=0:19;
 u2=0:9
 c2=c(0.5, 1.0, 1.5)
 c2*round(sqrt(0.1),2) # true sill=0.1
 
 ##---------------------------------------------------##
 ## True Probability of Agreement: Separable Case

 # Examples
 Psi_ST(h=0.1, u=1, C=0.5, paramcorr=param0.exp_exp, corrmodel="Exp_Exp", trend_t="linear")
 Psi_ST(h=0.1, u=1, C=0.5, paramcorr=param0.exp_exp[1,], corrmodel="Exp_Exp", trend_t="linear")
 Psi_ST(h=0.1, u=1, C=0.5, paramcorr=param0.exp_exp[2,], corrmodel="Exp_Exp", trend_t="linear")
 
 # Calculation for all cases
 Psi_ST_separable = list()
 for(k in 1: nrow(param0.exp_exp) ){
    aux = array(NA, dim=c(length(h2), length(u2), length(c2)) )
    for(i in 1:length(h2)){
       for(j in 1:length(u2)){
          for(l in 1:length(c2)){
             aux[i, j, l] <- Psi_ST(h=h2[i], u=u2[j], C=c2[l], paramcorr=param0.exp_exp[k,], corrmodel="Exp_Exp", trend_t="linear")
          }
       }
    }
    Psi_ST_separable[[k]] <- aux
 }
 dim(Psi_ST_separable[[1]]) # h:dim1, u:dim2, c:dim3

 ##-------------------------------------------##
 ## Estimated values of PA for NS=20
 Est_Psi_ST_NS20_T10_separable = list()
 for(k in 1:length(fit2_NS20_T10_separable)){
   aux = array(NA, dim=c(nrow(fit2_NS20_T10_separable[[k]]), length(h2), length(u2), length(c2)) )
   for(i in 1:length(h2)){
     for(j in 1:length(u2)){
       for(l in 1:length(c2)){
         aux[, i, j, l] <- Psi_est_st(h=h2[i], u=u2[j], C=c2[l], param.est=fit2_NS20_T10_separable[[k]], corrmodel="Exp_Exp")
       }
     }
   }
   Est_Psi_ST_NS20_T10_separable[[k]] <- aux
 }
 dim(Est_Psi_ST_NS20_T10_separable[[1]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4
 dim(Est_Psi_ST_NS20_T10_separable[[2]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4

 ##-------------------------------------------##
 ## Estimated values of PA for NS=50
 Est_Psi_ST_NS50_T10_separable = list()
 for(k in 1:length(fit2_NS50_T10_separable)){
    aux = array(NA, dim=c(nrow(fit2_NS50_T10_separable[[k]]), length(h2), length(u2), length(c2)) )
    for(i in 1:length(h2)){
       for(j in 1:length(u2)){
          for(l in 1:length(c2)){
             aux[, i, j, l] <- Psi_est_st(h=h2[i], u=u2[j], C=c2[l], param.est=fit2_NS50_T10_separable[[k]], corrmodel="Exp_Exp")
          }
       }
    }
    Est_Psi_ST_NS50_T10_separable[[k]] <- aux
 }
 dim(Est_Psi_ST_NS50_T10_separable[[1]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4
 dim(Est_Psi_ST_NS50_T10_separable[[2]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4
 
 ##------------------------------------------##
 ## boxplot of results
 
 pos.u=1:3; pos.c=1:3 # u=c(0,1,2); C=c(0.5,1.0,1.5) 
 pos.h=c(0,2,4,6,8)+1; # h=c(0,1,3,10); 
 
 png("psi_st_separable1_ns20.png", width=1152, height=1024, res=150 )
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
   for(l in 1:3){
     titulo=paste0("c=",c2[l], ", u=", u2[j])
     boxplot(Est_Psi_ST_NS20_T10_separable[[1]][,pos.h, j, l],
             ylim=c(min(Est_Psi_ST_NS20_T10_separable[[1]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
             names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
     points(Psi_ST_separable[[1]][pos.h, j, l], pch=20, col=2, cex=2) 
     points(apply(Est_Psi_ST_NS20_T10_separable[[1]][,pos.h, j, l],2,mean), pch=20, col=4, cex=2) 
     legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
   }
 }
 dev.off()
 
 png("psi_st_separable2_ns20.png", width=1152, height=1024, res=150 )
 pos.h=c(0,2,4,6,8)+1;
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
   for(l in 1:3){
     titulo=paste0("c=",c2[l], ", u=", u2[j])
     boxplot(Est_Psi_ST_NS20_T10_separable[[2]][,pos.h, j, l], 
             ylim=c(min(Est_Psi_ST_NS20_T10_separable[[2]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
             names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
     points(Psi_ST_separable[[2]][pos.h, j, l], pch=20, col=2, cex=2) 
     points(apply(Est_Psi_ST_NS20_T10_separable[[2]][,pos.h, j, l],2,mean), pch=20, col=4, cex=2) 
     legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
   }
 }
 dev.off()
 
 png("psi_st_separable1_ns50.png", width=1152, height=1024, res=150 )
 pos.h=c(0,2,4,6,8)+1;
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
    for(l in 1:3){
       titulo=paste0("c=",c2[l], ", u=", u2[j])
       boxplot(Est_Psi_ST_NS50_T10_separable[[1]][,pos.h, j, l], 
               ylim=c(min(Est_Psi_ST_NS50_T10_separable[[1]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
               names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
       points(Psi_ST_separable[[1]][pos.h, j, l], pch=20, col=2, cex=2) 
       points(apply(Est_Psi_ST_NS50_T10_separable[[1]][,pos.h, j, l],2,mean), pch=20, col=4, cex=2) 
       legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
    }
 }
 dev.off()
 
 png("psi_st_separable2_ns50.png", width=1152, height=1024, res=150 )
 pos.h=c(0,2,4,6,8)+1;
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
    for(l in 1:3){
       titulo=paste0("c=",c2[l], ", u=", u2[j])
       boxplot(Est_Psi_ST_NS50_T10_separable[[2]][,pos.h, j, l], 
               ylim=c(min(Est_Psi_ST_NS50_T10_separable[[2]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
               names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
       points(Psi_ST_separable[[2]][pos.h, j, l], pch=20, col=2, cex=2) 
       points(apply(Est_Psi_ST_NS50_T10_separable[[2]][,pos.h, j, l],2,mean), pch=20, col=4, cex=2) 
       legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
    }
 }
 dev.off()
 

 #######################################################
 ##---------------------------------------------------##
 ## True Probability of Agreement: Nonseparable Case
 
 # Examples
 Psi_ST(h=0, u=1, C=0.5, paramcorr=param0.iacocesare, corrmodel="Iacocesare", trend_t="linear")
 Psi_ST(h=0, u=1, C=0.5, paramcorr=param0.iacocesare[1,], corrmodel="Iacocesare", trend_t="linear")
 Psi_ST(h=0, u=1, C=0.5, paramcorr=param0.iacocesare[2,], corrmodel="Iacocesare", trend_t="linear")
 
 Psi_est_st(h=0, u=0, C=0.5, param.est=fit2_NS20_T10_nonsepar[[1]], corrmodel="Iacocesare")
 Psi_est_st(h=0, u=0, C=0.5, param.est=fit2_NS50_T10_nonsepar[[1]], corrmodel="Iacocesare")
 
 round(fit2_NS20_T10_nonsepar[[1]][ Psi_est_st(h=0, u=0, C=0.5, param.est=fit2_NS20_T10_nonsepar[[1]], corrmodel="Iacocesare")!=1,],6)
 round(fit2_NS20_T10_nonsepar[[2]][ Psi_est_st(h=0, u=0, C=0.5, param.est=fit2_NS20_T10_nonsepar[[2]], corrmodel="Iacocesare")!=1,],6)

 sigmaD_ST(h=0, u=0, corrmodel="Iacocesare", paramcorr=data.frame(fit2_NS20_T10_nonsepar[[1]][1,], nugget=0)) 
 sigmaD_ST(h=0, u=0, corrmodel="Iacocesare", paramcorr=data.frame(fit2_NS20_T10_nonsepar[[1]][8,], nugget=0))
 
 # Calculation for all cases
 Psi_ST_nonseparable = list()
 for(k in 1:nrow(param0.iacocesare)){
   aux = array(NA, dim=c(length(h2), length(u2), length(c2)) )
   for(i in 1:length(h2)){
     for(j in 1:length(u2)){
       for(l in 1:length(c2)){
         aux[i, j, l] <- Psi_ST(h=h2[i], u=u2[j], C=c2[l], paramcorr=param0.iacocesare[k,], corrmodel="Iacocesare", trend_t="linear")
       }
     }
   }
   Psi_ST_nonseparable[[k]] <- aux
 }
 dim(Psi_ST_nonseparable[[1]]) # h:dim1, u:dim2, c:dim3

 ##-------------------------------------------##
 ## Estimated values of PA for NS=20
 Est_Psi_ST_NS20_T10_nonseparable = list()
 for(k in 1:length(fit4_NS20_T10_nonsepar)){
   aux = array(NA, dim=c(nrow(fit4_NS20_T10_nonsepar[[k]]), length(h2), length(u2), length(c2)) )
     for(i in 1:length(h2)){
       for(j in 1:length(u2)){
         for(l in 1:length(c2)){
         aux[, i, j, l] <- Psi_est_st(h=h2[i], u=u2[j], C=c2[l], param.est=fit4_NS20_T10_nonsepar[[k]], corrmodel="Iacocesare")
       }
     }
   }
   Est_Psi_ST_NS20_T10_nonseparable[[k]] <- aux
 }
 dim(Est_Psi_ST_NS20_T10_nonseparable[[1]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4
 dim(Est_Psi_ST_NS20_T10_nonseparable[[2]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4

 ##------------------------------------------##
 ## Estimated values of PA for NS=50
 Est_Psi_ST_NS50_T10_nonseparable = list()
 for(k in 1:length(fit4_NS50_T10_nonsepar)){
    aux = array(NA, dim=c(nrow(fit4_NS50_T10_nonsepar[[k]]), length(h2), length(u2), length(c2)) )
    for(i in 1:length(h2)){
       for(j in 1:length(u2)){
          for(l in 1:length(c2)){
             aux[, i, j, l] <- Psi_est_st(h=h2[i], u=u2[j], C=c2[l], param.est=fit4_NS50_T10_nonsepar[[k]], corrmodel="Iacocesare")
          }
       }
    }
    Est_Psi_ST_NS50_T10_nonseparable[[k]] <- aux
 }
 dim(Est_Psi_ST_NS50_T10_nonseparable[[1]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4
 dim(Est_Psi_ST_NS50_T10_nonseparable[[2]]) # valid sim: dim1, h:dim2, u:dim3, c:dim4
 
 
 ##------------------------------------------##
 ## boxplot of results
 
 pos.u=1:3; pos.c=1:3 # u=c(0,1,2); C=c(0.5,1.0,1.5) 
 pos.h=c(0,2,4,6,8)+1; # h=c(0,1,3,10); 
 
 min(Est_Psi_ST_NS20_T10_nonseparable[[1]][,pos.h, j, l])*c(0.5,0.75,0.95)[l]
 
 png("psi_st_nonseparable1_ns20.png", width=1152, height=1024, res=150 )
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
   for(l in 1:3){
     titulo=paste0("c=",c2[l], ", u=", u2[j])
     boxplot(Est_Psi_ST_NS20_T10_nonseparable[[1]][,pos.h, j, l], 
             ylim=c(min(Est_Psi_ST_NS20_T10_nonseparable[[1]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
             names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
     points(Psi_ST_nonseparable[[1]][pos.h, j, l], pch=20, col=2, cex=2) 
     points(apply(Est_Psi_ST_NS20_T10_nonseparable[[1]][,pos.h, j, l],2, mean), pch=20, col=4, cex=2) 
     legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
   }
 }
 dev.off()
 
 png("psi_st_nonseparable2_ns20.png", width=1152, height=1024, res=150 )
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
   for(l in 1:3){
     titulo=paste0("c=",c2[l], ", u=", u2[j])
     boxplot(Est_Psi_ST_NS20_T10_nonseparable[[2]][,pos.h, j, l], 
             ylim=c(min(Est_Psi_ST_NS20_T10_nonseparable[[2]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
             names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
     points(Psi_ST_nonseparable[[2]][pos.h, j, l], pch=20, col=2, cex=2) 
     points(apply(Est_Psi_ST_NS20_T10_nonseparable[[2]][,pos.h, j, l],2,mean), pch=20, col=4, cex=2) 
     legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
   }
 }
 dev.off()
 
 png("psi_st_nonseparable1_ns50.png", width=1152, height=1024, res=150 )
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
   for(l in 1:3){
     titulo=paste0("c=",c2[l], ", u=", u2[j])
     boxplot(Est_Psi_ST_NS50_T10_nonseparable[[1]][,pos.h, j, l], 
             ylim=c(min(Est_Psi_ST_NS50_T10_nonseparable[[1]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
             names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
     points(Psi_ST_nonseparable[[1]][pos.h, j, l], pch=20, col=2, cex=2) 
     points(apply(Est_Psi_ST_NS50_T10_nonseparable[[1]][,pos.h, j, l],2,mean), pch=20, col=4, cex=2) 
     legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
   }
 }
 dev.off()
 
 png("psi_st_nonseparable2_ns50.png", width=1152, height=1024, res=150 )
 par(mfrow=c(3,3), mar=c(4,4.2,3,1))
 for(j in 1:3){
   for(l in 1:3){
     titulo=paste0("c=",c2[l], ", u=", u2[j])
     boxplot(Est_Psi_ST_NS50_T10_nonseparable[[2]][,pos.h, j, l], 
             ylim=c(min(Est_Psi_ST_NS50_T10_nonseparable[[2]][,pos.h, j, l])*c(0.7,0.96,0.999)[l],1),
             names=pos.h-1, xlab="||h||", ylab=expression(psi[c](h,u)), main=titulo, pch=20)
     points(Psi_ST_nonseparable[[2]][pos.h, j, l], pch=20, col=2, cex=2) 
     points(apply(Est_Psi_ST_NS50_T10_nonseparable[[2]][,pos.h, j, l],2,mean), pch=20, col=4, cex=2) 
     legend("bottomleft", legend=c("True", "Estimate"), col=c(2,4), pch=20, inset=0.01, bg="gray90")
   }
 }
 dev.off()