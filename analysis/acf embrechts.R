library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
load(paste0(datalib,"boot_embrechts5.Rda"))
load(paste0(datalib,"dowjones.Rda"))

pp   = dowjones$pp.neg
pars = dowjones$pars
Tmax = dowjones$Tmax

# These constraints ensures mu>0, n<1 and B>0:
# ui0 = diag(c(1,-1,1))
# ci0 = c(0,-1,0)
(res = mle(pp, Tmax, pars))

x = seq(0,Tmax,1)
L.exp = integrated_intensity(x, pp, res)
s.exp = interpol(cbind(x,L.exp), pp)

s = s.exp
w = diff(s)
v = w/mean(w);mean(w);sd(w)


layout(matrix(c(1,1,2,3,4,5), nr=3, nc=2, byrow=TRUE))
par(mar=c(3,3,1,1), oma=c(0,0,3,0))
plot_acf(diff(pp), lag=100)
mtext("ACF for Dow Jones")
# par(mfrow=c(2,2))

plot_acf(w, lag=100); mtext("Uncorrected waiting times", 3, line=0)
plot_acf(v, lag=100); mtext("Mean corrected waiting times", 3, line=0)
plot_acf(w^2, lag=100); mtext("Uncorrected waiting times squared", 3, line=0)
plot_acf(v^2, lag=100); mtext("Mean corrected waiting times squared", 3, line=0)
