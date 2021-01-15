library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
load(paste0(datalib,"boot_embrechts.Rda"))
load(paste0(datalib,"dowjones.Rda"))

pp   = dowjones$pp.neg
pars = dowjones$pars
Tmax = dowjones$Tmax

# These constraints ensures mu>0, n<1 and B>0:
ui0 = diag(c(1,-1,1))
ci0 = c(0,-1,0)
nrep = 100
pars
res = t(replicate(nrep,mle(pp, Tmax, c(runif(1,0,2*length(pp)/Tmax),
                                       runif(1,0,1),
                                       runif(1,0,2*length(pp)/Tmax)), dt=0, ui0=ui0,ci0=ci0)))

res = apply(res,2,mean)
log_lik(pp, Tmax, pars)
log_lik(pp, Tmax, res)

parnum=1
parname=expression(mu)
cols = add.alpha(c('steelblue', 'tomato2', 'tomato2'),.8)
cols = add.alpha(c('steelblue','tomato2'),.5)
plotpar <- function(parnum, parname="", cols=add.alpha('black',.5), yl=NULL){
  tmp = numeric(0)
  for(i in 1:4){
    tmp = cbind(tmp,boot[[i]][,parnum])  
  }
  if(is.null(yl)){
    yl = c(0,max(tmp))  
  }
  
  boxplot(tmp, col=cols[1], border='black', pch=16, cex=.75, xaxt='n', ylim=yl)
  abline(h=pars[parnum], lty=3, lwd=2)
  abline(h=res[parnum], lty=1, lwd=2)
  tmp = apply(tmp,2, make_ci)
  for(i in 1:4){
    segments(i,tmp[1,i],i,tmp[2,i], lwd=50, col=cols[2], lend='butt')
  }
  lbls = c("Parametric FIB", "Non-parametric FIB", "Parametric RIB", "Non-parametric RIB")
  text(1:4,.45*par("usr")[3], lbls, srt=0, xpd = NA, cex=.85)
  mtext(parname,side = 3, line=0)
}

fn = paste0(plotlib,"boot_embrechts.pdf")
pdf(fn, width=10, height=10)
    layout(1:3)
    par(mar=c(3,2,2,1), bty='n', oma=c(1,1,6,1))
    plotpar(1, expression(mu), cols, yl=c(0,.05))
    plotpar(2, "n", cols, yl=c(0,1))
    plotpar(3, expression(beta), cols, yl=c(0,.05))
    lgnd = c("MLE", "Sim. params", "Int-Q range","95% CI")
    legend("topleft",lgnd, col=c('black','black',cols), lty=c(1,3,NA,NA), pch=c(NA,NA,15,15), cex=1, bty='n', pt.cex = 2)
    mtext("Bootstrap Embrechts", 3, line=3, outer=TRUE, cex=1.5)  
    mtext("Using 399 bootstrap samples", 3, line=1.5, outer=TRUE)  
dev.off()
pars
res

