library(bootha)
library(misc)
library(plot3D)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
load(paste0(datalib,"dowjones.Rda"))

pp   = dowjones$pp.neg
pars = dowjones$pars
Tmax = dowjones$Tmax
res  = mle(pp, Tmax, pars,0)

ires = 50
res
mus  = seq(0.005,.025, length=ires)
ns   = seq(.5,1-1e-10,length=ires)
Bs   = seq(.005,.035, length=ires)
llks = array(dim=c(ires,ires,ires), dimnames = list(paste0("mu",1:ires),paste0("n",1:ires),paste0("B",1:ires)))
for(i in 1:ires){
  for(j in 1:ires){
    for(k in 1:ires){
      cat("\ri: ",ifelse(i<10," ",""),i, " j: ",ifelse(j<10," ",""),j," k: ",ifelse(k<10," ",""),k,"        ", sep="")  
      llks[i,j,k] = log_lik(pp, Tmax,c(mus[i],ns[j],Bs[k]))
    }
  }
}
# tmp  = llks
# llks = llks-min(llks)

idxm = which((mus-res[1])^2==min((mus-res[1])^2))
idxn = which((ns-res[2])^2==min((ns-res[2])^2))
idxB = which((Bs-res[3])^2==min((Bs-res[3])^2))

par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(0,2,0,2))
cols = colorRampPalette(c("dodgerblue","red"))(10000)
lbl = expression(bold("Likelihood function"~beta~"fixed"))
persp3D(mus,ns, llks[,,idxB], theta = -60, phi = 20, xlab="mu",ylab="n",col = cols, main="", colkey=FALSE, expand=.8, ticktype='detailed', facets = TRUE)
# points3D(res[1],res[2],log_lik(pp,Tmax,pars), add=TRUE, pch=16, col=add.alpha('black',.85))
scatter3D(res[1], res[2], log_lik(pp,Tmax,pars), add=TRUE, type='p', pch=124, col=add.alpha('black',.85), cex=2)
mtext(lbl, side=3, line=-10)

lbl = expression(bold("Likelihood function"~"n"~"fixed"))
persp3D(mus,Bs, llks[,idxn,], theta = 30, phi = 20, xlab="mu",ylab="B",col = cols, main="", colkey=FALSE, expand=.8, ticktype='detailed', facets = TRUE)
# points3D(res[1],res[3],log_lik(pp,Tmax,pars), add=TRUE, pch=16, col=add.alpha('black',.85), cex=2)
scatter3D(res[1], res[3], log_lik(pp,Tmax,pars), add=TRUE, type='p', pch=124, col=add.alpha('black',.85), cex=2)
mtext(lbl, side=3, line=-10)

lbl = expression(bold("Likelihood function"~mu~"fixed"))
persp3D(ns, Bs, llks[idxm,,], theta = -30, phi = 20, xlab="n",ylab="B",col = cols, main="", colkey=FALSE, expand=.8, ticktype='detailed', facets = TRUE)
# points3D(res[2],res[3],log_lik(pp,Tmax,pars), add=TRUE, pch=16, col=add.alpha('black',.85))
scatter3D(res[2], res[3], log_lik(pp,Tmax,pars), add=TRUE, type='p', pch=124, col=add.alpha('black',.85), cex=2)
mtext(lbl, side=3, line=-10)

round(res,3)
