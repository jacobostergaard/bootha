misc::clean_up()
library("PtProcess")
data("Phuket")
Phuket$magnitude <- Phuket$magnitude - 4.95 

dmagn_mark <- function(x, data, params) {
  if (params[7] > 0) {
    lambda <- etas_gif(data, x[, "time"], params = params[1:5])
    y <- dgamma(x[, "magnitude"], shape = 1 + sqrt(lambda) * params[7], rate = params[6], log = TRUE)
  } else y <- dexp(x[, "magnitude"], rate = params[6], log = TRUE) 
  
  return(y)

}

rmagn_mark <- function(ti, data, params) {
  if (params[7]>0) {
    lambda <- etas_gif(data, ti, params = params[1:5])
    y <- rgamma(1, shape = 1 + sqrt(lambda) * params[7], rate = params[6]) 
  } else y <- rexp(1, rate = params[6])
    
  return(list(magnitude = y))
}

TT <- c(0, 1827)
params <- c(0.05, 3.1, 1.3, 0.02, 1.1, 1/mean(Phuket$magnitude), 0) 

head(Phuket)
range(Phuket$time)

x <- mpp(data = Phuket, gif = etas_gif, mark = list(dmagn_mark, rmagn_mark), params = params, TT = TT, gmap = expression(params[1:5]), mmap = expression(params))


expmap <- function(y, p) {
    y$params[1:5] <- exp(p)
    return(y)
  }
initial <- log(params[1:5])
z <- optim(initial, neglogLik, object = x, pmap = expmap, control = list(trace = 1, maxit = 100))

initial <- z$par
z <- nlm(neglogLik, initial, object = x, pmap = expmap, print.level = 2, iterlim = 500, typsize = initial)
x0 <- expmap(x, z$estimate)


allmap <- function(y, p) {
    y$params <- exp(p)
    return(y)
  }
initial <- log(c(0.05, 3.1, 13, 0.02, 1.1, 1/mean(Phuket$magnitude), 0.1)) 

z <- optim(initial, neglogLik, object = x, pmap = allmap, control=list(trace = 1, maxit = 200)) 
initial <-z$par

z <- nlm(neglogLik, initial, object = x, pmap = allmap, print.level = 2, iterlim = 500, typsize = initial)
x1 <- allmap(x, z$estimate)



print(logLik(x0))
print(logLik(x1))






plot(residuals(x1), xlab = "Event Number", ylab = "Transformed Time", pty = "s")
points(residuals(x0), lty = 2, type = "l", col = "red")
big <- list(id = c(35, 511, 733, 875, 972, 1124, 1194), date = c("26Dec04", "28Mar05", "24Jul05", "16May06", "12Sep07", "20Feb08", "27Jun08"))
axis(3, at = big$id, labels = big$date, las = 2)
abline(v = big$id, lty = 3, col = "blue")
abline(a = 0, b = 1, lty = 3, col = "blue")

n <- nrow(Phuket)
plot(residuals(x1) - 1:n, xlab = "Event Number", ylab = "Cusum")
points(residuals(x0) - 1:n, lty = 2, type = "l", col = "red")
axis(3, at = big$id, labels = big$date, las = 2)
abline(v = big$id, lty = 3, col = "blue")
abline(a = 0, b = 0, lty = 3, col = "blue")


lambda <- etas_gif(data = x1$data, evalpts = x1$data$time, params = x1$params[1:5])
mean1 <- (1 + sqrt(lambda) * x1$params[7]) / x1$params[6]
plot(ts(cumsum(x1$data$magnitude - mean1)), ylab = "Cusum", xlab = "Event Number", ylim = c(-13, 20))

mean0 <- 1/x0$params[6]
points(ts(cumsum(x0$data$magnitude - mean0)), type = "l", lty = 2, col = "red")
abline(h = 0, lty = 3, col = "blue")
axis(3, at = big$id[-1], labels = big$date[-1], cex.axis = 0.7, las = 2)
abline(v = big$id[-1], lty = 3, col = "blue")






stop.cond <- function(data) {
    n <- nrow(data)
    return(data$magnitude[n] >= 1.55)
  }
x2 <- x1
x2$TT <- c(1827, Inf)
y <- rep(NA, 2000)


for (i in 1:2000) {
    print(i)
    x3 <- simulate(x2, seed = i, stop.cond = stop.cond)
    y[i] <- x3$data$time[nrow(x3$data)]
  }
hist(y - 1827, breaks = seq(0, max(y - 1827), length.out = 20), xlab = "Days Since 1 January 2009", main = "")
z <- quantile(y - 1827, probs = c(0.99, 0.95, 0.90, 0.80, 0.50)) 
abline(v = z, lty = 2, col = "blue")
axis(3, at = z, labels = c(0.99, 0.95, 0.90, 0.80, 0.50))
box()