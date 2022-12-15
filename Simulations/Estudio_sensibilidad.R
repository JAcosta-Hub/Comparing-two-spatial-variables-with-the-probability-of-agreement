###############################

 library(RandomFields)
 library(geoR)
 source("Functions_prob_agreement.R")

########################################
# variance of the difference for the Matern Covariance

 h=seq(0,1.5,by=0.01)
 M1=sigmaD(h, sigma2=c(1,1), phi=c(0.2,0.2,0.2), nu=c(0.5,0.5,0.5), rhoxy=0)
 M2=sigmaD(h, sigma2=c(1,1), phi=c(0.2,0.2,0.2), nu=c(0.5,0.5,0.5), rhoxy=0.25)
 M3=sigmaD(h, sigma2=c(1,1), phi=c(0.2,0.2,0.2), nu=c(0.5,0.5,0.5), rhoxy=0.5)
 M4=sigmaD(h, sigma2=c(1,1), phi=c(0.2,0.2,0.2), nu=c(0.5,0.5,0.5), rhoxy=0.75)
 M5=sigmaD(h, sigma2=c(1,1), phi=c(0.2,0.2,0.2), nu=c(0.5,0.5,0.5), rhoxy=1)

 plot(h, M1, type="l", ylim=c(0,2))
 points(h, M2, type="l", col=2)
 points(h, M3, type="l", col=3)
 points(h, M4, type="l", col=4)
 points(h, M5, type="l", col=5)

 
###################################################### 
## Probability of Agreement
 
 Mpsi.rho1=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0)
 Mpsi.rho2=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.25)
 Mpsi.rho3=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.rho4=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.75)
 Mpsi.rho5=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=1) 

 Mpsi.muD1=Psi(h=h, C=2, muD=0.0, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.muD2=Psi(h=h, C=2, muD=0.25, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.muD3=Psi(h=h, C=2, muD=0.50, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.muD4=Psi(h=h, C=2, muD=0.75, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.muD5=Psi(h=h, C=2, muD=1.00, sigma2=c(1,1), phi=c(0.2,0.2,0.2), rhoxy=0.5) 
 
 Mpsi.sigma1=Psi(h=h, C=2, muD=0, sigma2=c(1,0.8), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.sigma2=Psi(h=h, C=2, muD=0, sigma2=c(1,0.9), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.sigma3=Psi(h=h, C=2, muD=0, sigma2=c(1,1.0), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.sigma4=Psi(h=h, C=2, muD=0, sigma2=c(1,1.1), phi=c(0.2,0.2,0.2), rhoxy=0.5)
 Mpsi.sigma5=Psi(h=h, C=2, muD=0, sigma2=c(1,1.2), phi=c(0.2,0.2,0.2), rhoxy=0.5)

 Mpsi.phi1=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.10), rhoxy=0.5)
 Mpsi.phi2=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.15), rhoxy=0.5)
 Mpsi.phi3=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.20), rhoxy=0.5)
 Mpsi.phi4=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.25), rhoxy=0.5)
 Mpsi.phi5=Psi(h=h, C=2, muD=0, sigma2=c(1,1), phi=c(0.2,0.2,0.30), rhoxy=0.5)
 
###################################################### 
## Plotting Agreement Probabilities

 png("psi_matern_rho.png", width=768, height=512, res=130 )
 par(mar=c(4.5, 4.5, 1.5, 1.5))
 leg.rho=expression(rho[XY]==0.00, rho[XY]==0.25, rho[XY]==0.50, rho[XY]==0.75, rho[XY]==1.00)
 plot(h, Mpsi.rho1, type="l", ylim=c(0.8,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)) ) 
 points(h, Mpsi.rho2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.rho3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.rho4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.rho5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg.rho, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 dev.off()
 
 png("psi_matern_muD.png", width=768, height=512, res=130 )
 par(mar=c(4.5, 4.5, 1.5, 1.5))
 leg.muD=expression(mu[D]==0.00, mu[D]==0.25, mu[D]==0.50, mu[D]==0.75, mu[D]==1.00)
 plot(h, Mpsi.muD1, type="l", ylim=c(0.7,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)) ) 
 points(h, Mpsi.muD2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.muD3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.muD4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.muD5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg.muD, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 dev.off()
 
 
 png("psi_matern_sigma.png", width=768, height=512, res=130 )
 par(mar=c(4.5, 4.5, 1.5, 1.5))
 leg.sigY=expression(sigma[Y]==0.8, sigma[Y]==0.9, sigma[Y]==1.0, sigma[Y]==1.1, sigma[Y]==1.2)
 plot(h, Mpsi.sigma1, type="l", ylim=c(0.8,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)) ) 
 points(h, Mpsi.sigma2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.sigma3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.sigma4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.sigma5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg.sigY, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 dev.off()
 
 
 png("psi_matern_phi.png", width=768, height=512, res=130 )
 par(mar=c(4.5, 4.5, 1.5, 1.5))
 leg.phi=expression(phi[XY]==0.10, phi[XY]==0.15, phi[XY]==0.20, phi[XY]==0.25, phi[XY]==0.30)
 plot(h, Mpsi.phi1, type="l", ylim=c(0.8,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)) ) 
 points(h, Mpsi.phi2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.phi3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.phi4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.phi5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg.phi, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 dev.off()

  

###################################################### 
## Probability of Agreement
 
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho1=Psi(h=h, C=2, muD=0, rhoxy=0)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho2=Psi(h=h, C=2, muD=0, rhoxy=0.25)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho3=Psi(h=h, C=2, muD=0, rhoxy=0.5)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho4=Psi(h=h, C=2, muD=0, rhoxy=0.75)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho5=Psi(h=h, C=2, muD=0, rhoxy=1)

 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho1=Psi(h=h, C=2, muD=1, rhoxy=0)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho2=Psi(h=h, C=2, muD=1, rhoxy=0.25)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho3=Psi(h=h, C=2, muD=1, rhoxy=0.5)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho4=Psi(h=h, C=2, muD=1, rhoxy=0.75)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho5=Psi(h=h, C=2, muD=1, rhoxy=1)
 
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho1=Psi(h=h, C=2, muD=0, sigma2=c(1,2), rhoxy=0)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho2=Psi(h=h, C=2, muD=0, sigma2=c(1,2), rhoxy=0.25)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho3=Psi(h=h, C=2, muD=0, sigma2=c(1,2), rhoxy=0.5)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho4=Psi(h=h, C=2, muD=0, sigma2=c(1,2), rhoxy=0.75)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho5=Psi(h=h, C=2, muD=0, sigma2=c(1,2), rhoxy=1)
 
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho1=Psi(h=h, C=2, muD=1, sigma2=c(1,2), rhoxy=0)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho2=Psi(h=h, C=2, muD=1, sigma2=c(1,2), rhoxy=0.25)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho3=Psi(h=h, C=2, muD=1, sigma2=c(1,2), rhoxy=0.5)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho4=Psi(h=h, C=2, muD=1, sigma2=c(1,2), rhoxy=0.75)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho5=Psi(h=h, C=2, muD=1, sigma2=c(1,2), rhoxy=1)
 
 
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY2.rho1=Psi(h=h, C=2, muD=0, phi=c(0.2,0.2,0.5), rhoxy=0)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY2.rho2=Psi(h=h, C=2, muD=0, phi=c(0.2,0.2,0.5), rhoxy=0.25)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY2.rho3=Psi(h=h, C=2, muD=0, phi=c(0.2,0.2,0.5), rhoxy=0.5)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY2.rho4=Psi(h=h, C=2, muD=0, phi=c(0.2,0.2,0.5), rhoxy=0.75)
 Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY2.rho5=Psi(h=h, C=2, muD=0, phi=c(0.2,0.2,0.5), rhoxy=1)
 
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY2.rho1=Psi(h=h, C=2, muD=1, phi=c(0.2,0.2,0.5), rhoxy=0)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY2.rho2=Psi(h=h, C=2, muD=1, phi=c(0.2,0.2,0.5), rhoxy=0.25)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY2.rho3=Psi(h=h, C=2, muD=1, phi=c(0.2,0.2,0.5), rhoxy=0.5)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY2.rho4=Psi(h=h, C=2, muD=1, phi=c(0.2,0.2,0.5), rhoxy=0.75)
 Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY2.rho5=Psi(h=h, C=2, muD=1, phi=c(0.2,0.2,0.5), rhoxy=1)
 
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY2.rho1=Psi(h=h, C=2, muD=0, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY2.rho2=Psi(h=h, C=2, muD=0, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0.25)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY2.rho3=Psi(h=h, C=2, muD=0, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0.5)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY2.rho4=Psi(h=h, C=2, muD=0, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0.75)
 Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY2.rho5=Psi(h=h, C=2, muD=0, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=1)
 
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY2.rho1=Psi(h=h, C=2, muD=1, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY2.rho2=Psi(h=h, C=2, muD=1, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0.25)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY2.rho3=Psi(h=h, C=2, muD=1, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0.5)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY2.rho4=Psi(h=h, C=2, muD=1, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=0.75)
 Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY2.rho5=Psi(h=h, C=2, muD=1, sigma2=c(1,2), phi=c(0.2,0.2,0.5), rhoxy=1)
 
 
 
 
 leg=expression(rho[XY] == 0, rho[XY] == 0.25, rho[XY] == 0.5, rho[XY] == 
                  0.75, rho[XY] == 1)
 plot(h, Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho1, type="l",
      ylim=c(0.6,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)), 
      main=expression(mu[D]==0, sigma[X]==1, sigma[Y]==1) ) 
 points(h, Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.muD0.sigX1.sigY1.phiX1.phiY1.phiXY1.rho5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 
 plot(h, Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho1, type="l",
      ylim=c(0.6,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)), 
      main=expression(mu[D]==1, sigma[X]==1, sigma[Y]==1) ) 
 points(h, Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.muD1.sigX1.sigY1.phiX1.phiY1.phiXY1.rho5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 
 plot(h, Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho1, type="l",
      ylim=c(0.6,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)), 
      main=expression(mu[D]==0, sigma[X]==1, sigma[Y]==2) ) 
 points(h, Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.muD0.sigX1.sigY2.phiX1.phiY1.phiXY1.rho5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 
 
 plot(h, Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho1, type="l",
      ylim=c(0.6,1.02), lwd=1.5, xlab="||h||", ylab=expression(psi[c](h)), 
      main=expression(mu[D]==1, sigma[X]==1, sigma[Y]==2) ) 
 points(h, Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho2, type="l", col=2, lwd=1.5)
 points(h, Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho3, type="l", col=3, lwd=1.5)
 points(h, Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho4, type="l", col=4, lwd=1.5)
 points(h, Mpsi.muD1.sigX1.sigY2.phiX1.phiY1.phiXY1.rho5, type="l", col=5, lwd=1.5)
 legend("topright", legend=leg, lty=1, lwd=2, col=1:5, inset = 0.02, bg="lightgray")
 

 