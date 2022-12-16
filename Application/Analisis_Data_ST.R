##########################################
##--------------------------------------##
## Librerias y fuentes

 library(GeoModels)
 library(geoR)
 source("Functions_prob_agreement.R")
 source("Adaptation_plot_variograms_st.R")
 load("Imagenes.RData")
 datos_gcc <- read.csv("image_Gcc_ROI.csv", header=TRUE, skip = 17)

 ##-------------------------------------##
 ## G_cc calculation
 images_st <- gcc_st(raster15_images)

 plot(images_st$coords)

 ##-------------------------------------##
 ## Empirical Variogram spacetime
 
 vario_st = GeoVariogram(data=images_st$data,
                         coordx=images_st$coords, 
                         coordt=images_st$times,
                         maxtime=10, maxdist=30)
 plot(vario_st)
 
 # Remove linear time trend
 reg1=lm(c(t(images_st$data))~0+images_st$X)
 res1=matrix(reg1$residuals, length(images_st$times), 
             nrow(raster15_images[[1]])*ncol(raster15_images[[1]]), 
             byrow=TRUE)
 
 plot(c(t(images_st$data))~images_st$X[,2], xlab="time", ylab="gcc",
      main=paste0("slope = ",round(reg1$coefficients[2],4), " & ", "intercept = ", round(reg1$coefficients[1],3)))
 abline(coef=reg1$coefficients, col=2, lwd=2)

 vario1_st = GeoVariogram(data=res1, 
                          coordx=images2_st$coords, 
                          coordt=images2_st$times,
                          maxtime=10, maxdist=30)
 
 png("variograms_st.png", width=1024, height=1024, res=150)
 plot.GeoVariogram2(vario1_st)
 dev.off()

 ############################################################
 ##--------------------------------------------------------##
 ## Estimations

 ##--------------------------------------------------------##
 ## Exponential separable model
 
 ## Model, Initial and fixed parameters
 corrmodel="exp_exp"
 CorrParam("exp_exp")
 
 fixed=list(nugget=0)
 start=list(mean=0.5, mean1=0.01,
            sill=0.001, scale_s=20/log(1/0.05), scale_t=1)
 
 ## Fitting 
 fit1.exp_exp <- GeoFit(data=images_st$data,
                           coordx=images_st$coords, 
                           coordt=images_st$times,
                           X=images_st$X,
                           corrmodel=corrmodel, 
                           maxtime=2, maxdist=start$scale_s/2,
                           likelihood="Marginal", type="Pairwise",
                           start=start, fixed=fixed)
 fit1.exp_exp

 
 ##--------------------------------------------------------##
 ## Iacocesare nonseparable model
 
 ## Model, Initial and fixed parameters
 corrmodel="Iacocesare"
 CorrParam("Iacocesare")
 
 fixed=list(nugget=0)
 start=list(mean=0.5, mean1=0.01,
            sill=0.001, scale_s=20/log(1/0.05), scale_t=1,
            power2=2, power_s=1, power_t=1)
 
 ## Fitting 
 fit1.iacocesare <- GeoFit(data=images2_st$data,
                           coordx=images2_st$coords, 
                           coordt=images2_st$times,
                           X=images2_st$X,
                           corrmodel=corrmodel, 
                           maxtime=2, maxdist=start$scale_s/2,
                           likelihood="Marginal", type="Pairwise",
                           start=start, fixed=fixed)
 fit1.iacocesare
 
 ##------------------------------##
 ## Estimated trend graphs
 
 png("GCCs.png", width = 1024, height = 768, res=130)
 par(mar=c(4.2,4.5,2,1))
 plot(c(t(images_st$data))~images_st$X[,2], xlab="time", ylab=expression(G[cc]), 
      ylim=c(0.0, 0.6), pch=20, col="gray70")
 points(1:15, datos_gcc$gcc, pch=20, cex=1.75)
 legend("topleft", legend=expression(G[cc]), inset = 0.05, pch=20, pt.cex=1.75)
 legend("bottomright", lwd=2, lty=1, col=2:4, inset = 0.05,  bg="gray80",
        legend=c(
           paste0("Linear model, slope = ",round(reg1$coefficients[2],5)," & ", 
                  "intercept = ",round(reg1$coefficients[1],4)),
           paste0("Exponential, slope = ",round( fit1.exp_exp$param$mean1,5)," & ",
                  "intercept = ",round(fit1.exp_exp$param$mean,4)),
           paste0("Iacocesare, slope = ",round( fit1.iacocesare$param$mean1,5)," & ", 
                  "intercept = ",round(fit1.iacocesare$param$mean,4)) )        )
 abline(coef=reg1$coefficients, col=2, lwd=2)
 abline(coef= unlist(fit1.exp_exp$param)[1:2], col=3, lwd=2)
 abline(coef= unlist(fit1.iacocesare$param)[1:2], col=4, lwd=2)
 dev.off()

 ##---------------------------------------##
 ## Model selection
 
 print(booktabs=TRUE,  xtable::xtable(
    data.frame(
       logCompLik=rbind(exp_exp=fit1.exp_exp$logCompLik,
                        iacocesare=fit1.iacocesare$logCompLik),
       AIC=rbind(
          exp_exp=-2*fit1.exp_exp$logCompLik+2*length(fit1.exp_exp$param),
          iacocesare=-2*fit1.iacocesare$logCompLik+2*length(fit1.iacocesare$param)))))
 
 
 which.max( c(
    fit1.exp_exp$logCompLik,
    fit1.iacocesare$logCompLik))
 
 which.min( c(
    -2*fit1.exp_exp$logCompLik+2*length(fit1.exp_exp$param),
    -2*fit1.iacocesare$logCompLik+2*length(fit1.iacocesare$param)) )
 # Iacocesare is selected 

 ##------------------------------##
 ## Result of model selected
 
 ## Parameters Estimation
 fixed=list(nugget=0)
 start=list(mean=0.5, mean1=0.01,
            sill=0.001, scale_s=20/log(1/0.05), scale_t=1,
            power2=2, power_s=1, power_t=1)
 
 print(booktabs=TRUE, xtable::xtable(
    rbind(
       initial=unlist(start[names(fit1.iacocesare$param)]),
       estimated=unlist(fit1.iacocesare$param)), digits=5))
 
 ## Practical Range
 u3=seq(0,10, by=0.1)
 h3=fit1.iacocesare$param$scale_s*(20^(1/fit1.iacocesare$param$power2)-1-(u3/fit1.iacocesare$param$scale_t)^fit1.iacocesare$param$power_t )^(1/fit1.iacocesare$param$power2)
 plot(u3, h3, type="l" )
 h3[1] 
 

 ##--------------------------------##
 ## Probabilidad de acuerdo
 
 h2=seq(0, 10, by=0.5);                      # h variation
 u2=0:5                                      # u variation
 c2=round(c(0.5, 1, 1.5, 2)*sqrt(0.00044),2) # c variation
 

 ##---------------------------------------------------##
 ## Probability of Agreement: Separable Case

 Psi_ST_separable = array(NA, dim=c(length(h2), length(u2), length(c2)) )
    for(i in 1:length(h2)){
       for(j in 1:length(u2)){
          for(l in 1:length(c2)){
             Psi_ST_separable[i, j, l] <- Psi_est_st(h=h2[i], u=u2[j], C=c2[l], 
                                                     param.est=as.data.frame((fit1.exp_exp$param)), corrmodel="Exp_Exp")
          }
       }
    }
 # Plot example
 j=1
 plot(h2, Psi_ST_separable[,j,1], type="l", col=1, main=paste0("u=",u2[j]), ylim=c(0,1))
 points(h2, Psi_ST_separable[,j,2], type="l", col=2)
 points(h2, Psi_ST_separable[,j,3], type="l", col=3)
 points(h2, Psi_ST_separable[,j,4], type="l", col=4)
 
 ##---------------------------------------------------##
 ## Probability of Agreement: Nonseparable Case

 Psi_ST_nonseparable = array(NA, dim=c(length(h2), length(u2), length(c2)) )
 for(i in 1:length(h2)){
    for(j in 1:length(u2)){
       for(l in 1:length(c2)){
          Psi_ST_nonseparable[i, j, l] <- Psi_est_st(h=h2[i], u=u2[j], C=c2[l], 
                                                     param.est=as.data.frame((fit1.iacocesare$param)),
                                                     corrmodel="Iacocesare")
       }
    }
 }
 # A 3-dimensional array (h, u and c), i.e, 
 # row is h variation, column is time variation, for each c value.

 png("PA_st.png", width=1024, height=768, res=130)
 par(mfrow=c(2,2), mar=c(4.5,4.5,2.5,1))
 for(j in 1:4){
    plot(h2, Psi_ST_nonseparable[,j,1], type="l", col=1, main=paste0("Time difference, u=",u2[j]), 
         xlim=c(0,16), ylim=c(0,1), xlab="||h||", ylab=expression(Psi[c](h,u)))
    points(h2, Psi_ST_nonseparable[,j,2], type="l", col=2)
    points(h2, Psi_ST_nonseparable[,j,3], type="l", col=3)
    points(h2, Psi_ST_nonseparable[,j,4], type="l", col=4)
    legend("bottomright", lwd=2, lty=1, col=1:4, inset = 0.05,  bg="gray80", 
           legend=c(paste0("c=",c2[1]), 
                    paste0("c=",c2[2]), 
                    paste0("c=",c2[3]), 
                    paste0("c=",c2[4]) )  )
 }
 dev.off()
 
 
 # save the principal results
 save(images_st, vario_st, vario1_st, reg1, res1,  fit1.iacocesare,  fit1.exp_exp,
      Psi_ST_separable, Psi_ST_nonseparable, file="PA_ST_images.RData")
 