library("spam64")
library('LatticeKrig')
# ncobj<-ncdf4::nc_open('/glade/work/mleduc/glow output/jan30storm/fx2000_ne120pg3L273.001_1min.cam.h1.0001-01-30-00000_glowout.001.nc') #creating nc object
# oxy<-ncdf4::ncvar_get( ncobj,varid='ETA', start = c(2,2,56,13), count=c(-1,-1,5,1) )
# lev<-ncdf4::ncvar_get( ncobj,varid='lev',start = 56, count = 5)
# lat<-ncdf4::ncvar_get(ncobj,varid='lat', start = 2, count=-1)
# lon<-ncdf4::ncvar_get(ncobj,varid='lon', start = 2, count=-1)
inputs <- R.matlab::readMat('data for lkrig.mat')
data <- inputs$data - rowMeans(inputs$data)
marginalStdevs <- matrixStats::rowSds(inputs$data)
ndata <- data/marginalStdevs 
gridpts <- cbind( inputs$LON-180,inputs$LAT )
############# first level
LKinfo1 <- LatticeKrig::LKrigSetup(gridpts ,LKGeometry="LKSphere",startingLevel=1,nlevel=3,a.wght=1.01,lambda=0.0001,alpha=c(1,0.5,0.25)/1.75)
fitMLElambda <- LKrigFindLambda(gridpts,ndata,LKinfo=LKinfo1)
LKinfo1 <-LKinfoUpdate(LKinfo1,lambda=fitMLElambda$summary[[6]])

# aWghtGrid <- 1+1/c(seq(0.1,0.9,by=0.1)^2,seq( 1,101,by=4 )^2)
# llambda <-seq(-12,0,by=0.5)
# par.grid <- list(a.wght=aWghtGrid)
# fitMLE1<- LKrig.MLE(gridpts,ndata,LKinfo=LKinfo1,
#                      par.grid= par.grid )
fitMLEawght <- LKrigFindLambdaAwght(gridpts,ndata,LKinfo=LKinfo1)
LKinfo1 <-LKinfoUpdate(LKinfo1,lambda=fitMLEawght$summary[[6]],a.wght = fitMLEawght$summary[[7]],rho=fitMLEawght$summary[[5]])
# ndx <- which.min(fitMLE1$summary[,2])
Phi <- LKrig.basis(gridpts,LKinfo = LKinfo1)
Q <- LKrig.precision(LKinfo1)
# simulation <- LKrig.sim(gridpts,LKinfo=LKinfo1,M=100)
############ Calculate likelihoods
# ndcs[1,] <- 
# ndcs[2,] <- sample(1:64800,24800,replace=FALSE)

lambda<-fitMLEawght$summary[[6]]
sigmasq<-fitMLEawght$summary[[4]]^2
nll <- vector(length=5)
rho<-fitMLEawght$summary[[5]]
Q<- Q/rho
for(ii in 1:5){
  thesendcs <- sample(1:64800,24800,replace=FALSE)
  
  QpPtDiP <- Q+1/sigmasq*t(Phi[thesendcs,])%*%Phi[thesendcs,]
  Lq <- t(chol(spam2full(QpPtDiP)))
  projdata <- 1/sigmasq*t(Phi[thesendcs,])%*%ndata[thesendcs,]
  trbig<-sum((solve(Lq)%*%projdata)^2)/25
  trS <- 1/25*sum((ndata[thesendcs,]^2/(sigmasq)))
  ldd <- length(thesendcs)*log(sigmasq)
  ldq <- log(det(Q))
  ldqp <- 2*sum(log(diag(Lq)))
  nll[ii] <- ldqp + ldd + trS - trbig - ldq 
}




