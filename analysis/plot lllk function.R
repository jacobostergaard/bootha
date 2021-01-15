library(emhawkes)
library(hawkes)
library(ppboot)
library(misc)
library(parallel)
clean_up()

set.seed(1234)


mu   = 1
a    = 4
b    = 5
Tmax = 1000

set.seed(1234)
obs = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))


A = seq(3.5,6.5,.025)
B = seq(4.5,7.5,.025)
X = as.matrix(expand.grid(A,B))
X = lapply(seq_len(nrow(X)), function(i) X[i,])
system.time({
    Y = mclapply(X, function(x) if(x[1]<x[2]) llk(c(mu,x[1],x[2]), obs,package = "emhawkes"), mc.cores = 4)
    Y = unlist(lapply(Y, function(x) ifelse(is.null(x),NA,x )))
})

Z = cbind(as.matrix(expand.grid(A,B)),Y)
Z = matrix(Y,nr=length(A),nc=length(B))

mle = quiet(emhawkes::hfit(new("hspec", mu=mu, alpha=a, beta=b), diff(c(0,obs))))
sink()


deviations = seq(-.5,.5,.01)

llks = data.frame(mu=numeric(length(deviations)),a=numeric(length(deviations)),b=numeric(length(deviations)))

Z2 = Z
Z2[is.na(Z)]=-1
Z2[Z>0] = NA

Z1 = exp(Z-max(Z,na.rm = TRUE))



par(mfrow=c(1,2), bty='n', las=1, cex.lab=1.2)
lbl = expression(bold("Likelihood function"~mu~"fixed"))
image(B,A,t(Z1), xlab=expression(beta),ylab=expression(alpha),col = colorRampPalette(c("white","red"))(10000), main=lbl)
par(new=TRUE)
image(B,A,t(Z2), xlab=expression(beta),ylab=expression(alpha),col ='gray90')

abline(0,a/b, lty=3)
points(b,a, bg='dodgerblue', lwd=1, cex=2, pch=21)
points(mle$estimate[3],mle$estimate[2], lwd=1, bg='orange', cex=2, pch=21)
legend("bottomright", c("True","MLE"), col=c("dodgerblue","orange"), pch=16, cex=1, bty='n', x.intersp = .1, y.intersp = 3)

# lines(A-.05,B-.05, type='s')


for(i in 1:length(deviations)){
  tmp = deviations[i]
  llks$mu[i] = llk(c(mu*(1+tmp),a,b), obs)  
  llks$a[i] = llk(c(mu,a*(1+tmp),b), obs)  
  llks$b[i] = llk(c(mu,a,b*(1+tmp)), obs)  
}
lbl = "log-likelihoods with 1 varying parameter"
plot(0,0,type='n', xlim=range(100*deviations), ylim=range(llks)*1.1, xlab="% deviation from true value", ylab="log-likelihood", main=lbl)
lines(100*deviations, llks$mu, col='red', lwd=2)
lines(100*deviations, llks$a, col='dodgerblue', lwd=2)
lines(100*deviations, llks$b, col='black', lwd=2)
# abline(v=0+seq(-50,50,10), lty=3)
abline(h=llk(mle$estimate, obs), lty=2)
legend("topleft", c(expression(mu),expression(alpha),expression(beta), "Max log-likelihood"), col=c('red',"dodgerblue",'black', 'black'), lty=c(1,1,1,3), lwd=2, x.intersp = .1, seg.len = .4, bty='n', y.intersp = 3)


lbl = bquote(atop((mu~","~alpha~","~beta)=="("~.(mu)~","~.(a)~","~.(b)~")","T=[0,"~.(Tmax)~"]"))

mtext(lbl, side = 3, outer=TRUE, line=-3, cex=1.25)



a*1.05
b*1.05
mle$estimate
# mle$estimate/c(mu,a,b)
theta = c(mu, a, b)

c(0.52, 4.13, 4.28)/theta
c(1.06, 6.25, 7.20)/theta
c(0.98, 5.61, 6.85)/theta
c(1.00, 4.25, 5.15)/theta

a/b
4.13/4.28
6.25/7.20
5.61/6.85
4.25/5.15
