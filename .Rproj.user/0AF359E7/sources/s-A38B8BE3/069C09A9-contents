# Analysis of DJI
library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"

savePDF = FALSE
savePDF = TRUE

# Make point process of Dow Jones Index log return exceedances
    fn  = paste0(datalib,"DJI.csv")
    dji = read.csv(fn)
    dji$Date = as.Date(dji$Date)
    dji$log.return = c(NA,diff(log(dji$Close)))
    dji = dji[,c(1,5,7,8)] # relevant columns
    
    date.from  = as.Date("1994-01-01")
    date.to    = as.Date("2010-12-31")
    idx        = dji$Date>date.from & dji$Date < date.to # same period as in Embrechts paper
    names(dji) = c("date", "close", "volume", "log.return")
    dji        = dji[idx,]
    
    qts = quantile(dji$log.return, probs = c(.1,.9), na.rm = TRUE) 
    idx = which(dji$log.return < qts[1] | dji$log.return > qts[2]) # find days where the log return exceeds the .1 and .9 quantiles
    
    x    = as.numeric(dji$date[idx])
    x    = x-min(x)+1
    y    = dji$log.return[idx]
    
    pars = c(0.018,0.74,0.021) # Parameters from the Embrechts paper
    Tmax = max(x)
    dji  = list(pp.neg = as.numeric(x[y<0]), pp.pos = as.numeric(x[y>0]))
     
# Estimate MLE and rescaled time transformed waiting times
    pp   = dji$pp.neg # choose events of negative exceedances
    
    t_norm = 14
    t_norm = max(x)/length(pp)
  
    pp = pp/t_norm
    Tmax = Tmax/t_norm
    pars[c(1,3)]=pars[c(1,3)]*t_norm
  
    
    mle  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)

    
    dt = 0.1 # time resolution for integrating the intensity
    x  = seq(0,Tmax,dt)
    L  = integrated_intensity(x, pp, mle)
    s  = interpol(cbind(x,L), x=pp) # time transformed event times
    w  = c(s[1],diff(s))            # time transformed waiting times (first waiting time equals s[1])
    v  = w/mean(w);mean(w);sd(w)    # rescaled (mean corrected) time transformed waiting times
    
# Model diagnostics: evaluate model fit based on w and v
  
  fn  = paste0(plotlib,"dji_diagnostics.pdf")
  if(savePDF) pdf(file = fn, height=12, width=8)
    par(mfrow=c(3,2), mar=c(5,4,1,1), oma=c(0,0,0,0), bty='n')
    # Plot fitted model
        hi = Tmax+1
        lo = 0
        dt = (365.25*3/12)/t_norm  # time resolution to make a histogram like representation of the data (observed discrete time bins correspond to 1 day)
        x  = seq(0,Tmax,dt)
        y  = intensity(x, pp, mle)
        z  = pp[pp > lo & pp < hi]
        h  = hist(pp, breaks=seq(0,max(z)+dt,dt), plot=FALSE)
        
        # xlbl = paste("Days since", date.from)
        xlbl = "Year"
        ylbl = "Events per quarter (3 months)"
        col1 = add.alpha('steelblue',.5)
        col2 = add.alpha('tomato2',.75)
        col3 = add.alpha('black',.5)
        plot(h$breaks[-1], h$counts, type='l', xlim=c(lo,hi), bty='n', col=col1, lwd=3, xlab="", ylab="", xaxt='n')
        mtext(xlbl, side=1, line=2.5, cex=.75)
        mtext(ylbl, side=2, line=2.5, cex=.75)
        tmp = as.Date(paste0(seq(1995,2010,5),"-01-01"))
        xat = (as.numeric(tmp)-as.numeric(as.Date(date.from)))/t_norm
        axis(side = 1, at=xat, labels=substr(tmp,1,4))
        lines(x, y*dt, col=col2, lwd=2)
        points(z, rep(0,length(z)), pch=124, col=col3)
      
        legend("topleft", c("Observed intensity", "Fitted intensity", "Observed events"), pch=c(NA, NA, 124), lty=c(1,1,NA), col=c(col1, col2, col3), bty='n')
        
    # Plot s vs Poisson process bounds, this is in the Embrechts paper, but it is not too usefull...
        a = qnorm(.975)
        mu.tmp = length(s)/max(s)
        mu.hi  = mu.tmp + a*sqrt(mu.tmp/length(s))
        mu.lo  = mu.tmp - a*sqrt(mu.tmp/length(s))
        col1 = add.alpha('black',.75)
        col2 = add.alpha('tomato2',.75)
        col3 = add.alpha('steelblue',.5)
        xlbl = expression("Rescaled time"~Lambda(t))
        ylbl = expression("Time rescaled process"~ N(Lambda(t)))
        plot(s,1:length(s), type='n', xlab="", ylab="")
        mtext(xlbl, side=1, line=2.5, cex=.75)
        mtext(ylbl, side=2, line=2.5, cex=.75)
        yhi = (mu.hi-mu.tmp)*max(s)
        ylo = (mu.lo-mu.tmp)*max(s)
        Nt  = length(pp)
        polygon(c(0,Nt,Nt,0), c(yhi,Nt+yhi,Nt+ylo,ylo), border=NA, col=col3)
        abline(0,mu.tmp, col=col2)
        lines(s,1:length(s), type='s', col=col1)
        # abline((mu.lo-mu.tmp)*max(s),mu.tmp, lty=3)
        # abline((mu.lo-mu.tmp)*max(s),mu.tmp, lty=3)

      # Plot quantile-quantile plots of w and v
        xlbl = expression("Time transformed waiting times"~w[k]==Lambda(t[k]-t[k-1]))
        ylbl = "Theoretical"
        mlbl = "QQ plot"
        qqplot(w, rexp(10000, 1), xlim=c(0,8), ylim=c(0,8), col=add.alpha('red',.85), pch=16, ann=FALSE); abline(0,1)
        mtext(xlbl, side=1, line=2.5, cex=.75)
        mtext(ylbl, side=2, line=2.5, cex=.75)
        mtext(mlbl, side=3, line=0, cex=.75)
        
        xlbl = expression("Rescaled time transformed waiting times"~v[k]==w[k]/bar(w))
        qqplot(v, rexp(10000, 1), xlim=c(0,8), ylim=c(0,8), col=add.alpha('red',.85), pch=16, ann=FALSE); abline(0,1)
        mtext(xlbl, side=1, line=2.5, cex=.75)
        mtext(ylbl, side=2, line=2.5, cex=.75)
        mtext(mlbl, side=3, line=0, cex=.75)
        
      # Kolmogorov-Smirnov test and plots for the time transformed waiting times
        mlbl = expression("Kolmogorov-Smirnov plot of"~w)
        ks_plot(w)
        mtext(mlbl, side=3, line=0, cex=.75)
        mlbl = expression("Kolmogorov-Smirnov plot of"~v)
        ks_plot(v)
        mtext(mlbl, side=3, line=0, cex=.75)
  if(savePDF) dev.off()
        # There seems to be no worries so far.
        
# Plot autocorrelations
  fn  = paste0(plotlib,"dji_acf.pdf")
  if(savePDF) pdf(file = fn, height=12, width=8)
      layout(matrix(c(1,1,2,3,4,5), nr=3, nc=2, byrow=TRUE))
      par(mar=c(5,4,1,1), oma=c(0,0,1,0), bty='n')
      xlbl = "lag"
      ylbl = "ACF"
      mlbl = "Observed Dow Jones waiting times"
      plot_acf(diff(pp), lag=100, ann=FALSE)
      mtext(mlbl, side = 3, line=0, cex=.75)
      mtext(xlbl, side=1, line=2.5, cex=.75)
      mtext(ylbl, side=2, line=2.5, cex=.75)
      
      mlbl = expression("Time transformed waiting times"~w[k]==Lambda(t[k]-t[k-1]))
      plot_acf(w, lag=100, ann=FALSE)
      mtext(mlbl, side = 3, line=0, cex=.75)
      mtext(xlbl, side=1, line=2.5, cex=.75)
      mtext(ylbl, side=2, line=2.5, cex=.75)
      
      plot_acf(v, lag=100, ann=FALSE)
      mlbl = expression("Rescaled time transformed waiting times"~v[k]==w[k]/bar(w))
      mtext(mlbl, side = 3, line=0, cex=.75)
      mtext(xlbl, side=1, line=2.5, cex=.75)
      mtext(ylbl, side=2, line=2.5, cex=.75)
      
      plot_acf(w^2, lag=100, ann=FALSE)
      mlbl = expression("Time transformed waiting times"~w^2)
      mtext(mlbl, side = 3, line=0, cex=.75)
      mtext(xlbl, side=1, line=2.5, cex=.75)
      mtext(ylbl, side=2, line=2.5, cex=.75)
      
      plot_acf(v^2, lag=100, ann=FALSE)
      mlbl = expression("Rescaled time transformed waiting times"~v^2)
      mtext(mlbl, side = 3, line=0, cex=.75)
      mtext(xlbl, side=1, line=2.5, cex=.75)
      mtext(ylbl, side=2, line=2.5, cex=.75)
  if(savePDF) dev.off()
      # ACFs are ok, we assume independence of time transformed waiting times.
  
# Bootstrap samples
    cat(paste0("\nBegin at: ", Sys.time()),"\n")
      fn = paste0(datalib,"dji_boot.Rda")
      B = 399    # Set number of bootstrap samples
      ncores = 8 # Number of cpu cores to use in bootstrapping
      dji.boot = list() 
      system.time({
        dji.boot[[1]] = bootstrap_hawkes(pp, mle, Tmax, B, type='fib', parametric=TRUE, ncores=ncores)
      })
      system.time({
        dji.boot[[2]] = bootstrap_hawkes(pp, mle, Tmax, B, type='fib', parametric=FALSE, ncores=ncores)
      })
      system.time({
        dji.boot[[3]] = bootstrap_hawkes(pp, mle, Tmax, B, type='rib', parametric=TRUE, ncores=ncores)
      })
      system.time({
        dji.boot[[4]] = bootstrap_hawkes(pp, mle, Tmax, B, type='rib', parametric=FALSE, ncores=ncores)
      })
      names(dji.boot) = c("fib_par","fib_npar", "rib_par","rib_npar")
    cat(paste0("\nDone at: ", Sys.time()))
    
  
    
# Kolmogorov-Smirnov test    
    cat("\nKS test for time transformed waiting times w:\n")
    ks.test(w,'pexp')
    cat("\nKS test for rescaled (mean corrected) time transformed waiting times v:\n")
    ks.test(v,'pexp')
    
# Bootstrap Confidence Intervals
    CI = list()
    # Asypmtotic CIs
      I   = -solve(exp_hessian(pp, Tmax, mle))
      ses = sqrt(diag(I)) 
      CI$asymptotic = data.frame(mle=mle, ci.lo = mle-1.96*ses, ci.hi = mle+1.96*ses)
      
    # Three kinds of confidence intervals for all 4 types of bootstrap
      CI$naive$par_fib  = boot.ci(dji.boot$fib_par, pp, Tmax, pars, type  = 'naive')
      CI$naive$npar_fib = boot.ci(dji.boot$fib_npar, pp, Tmax, pars, type = 'naive')
      CI$naive$par_rib  = boot.ci(dji.boot$rib_par, pp, Tmax, pars, type  = 'naive')
      CI$naive$npar_rib = boot.ci(dji.boot$rib_npar, pp, Tmax, pars, type = 'naive')
      
      CI$tratios$par_fib  = boot.ci(dji.boot$fib_par, pp, Tmax, pars, type  = 'tratios')
      CI$tratios$npar_fib = boot.ci(dji.boot$fib_npar, pp, Tmax, pars, type = 'tratios')
      CI$tratios$par_rib  = boot.ci(dji.boot$rib_par, pp, Tmax, pars, type  = 'tratios')
      CI$tratios$npar_rib = boot.ci(dji.boot$rib_npar, pp, Tmax, pars, type = 'tratios')
      
      CI$bootstd$par_fib  = boot.ci(dji.boot$fib_par, pp, Tmax, pars, type  = 'bootstd')
      CI$bootstd$npar_fib = boot.ci(dji.boot$fib_npar, pp, Tmax, pars, type = 'bootstd')
      CI$bootstd$par_rib  = boot.ci(dji.boot$rib_par, pp, Tmax, pars, type  = 'bootstd')
      CI$bootstd$npar_rib = boot.ci(dji.boot$rib_npar, pp, Tmax, pars, type = 'bootstd')
      
    
      # xtable::xtable(CI$asymptotic, digits=3)
      # lapply(CI$naive, function(x) xtable::xtable(x, digits=3 ) )
      # lapply(CI$tratios, function(x) xtable::xtable(x, digits=3 ) )
      # lapply(CI$bootstd, function(x) xtable::xtable(x, digits=3 ) )
      
      round(CI$asymptotic,3)
      lapply(CI$naive, function(x) round(x,3))
      lapply(CI$tratios, function(x) round(x,3))
      lapply(CI$bootstd, function(x) round(x,3))
      # doption = options()$digits
      # options(digits=3)
      # options(digits=doption)