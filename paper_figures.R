# Plots for different dose-response models (Figure 1)
library(dreamer)
library(targets)
doses <- seq(from = 0, to = 10, length.out = 100)
d1 <- data.frame("dose"=doses)
# sigmoid EMAX
dat <- dreamer_data_emax(
  n_cohorts = rep(1, length(doses)),
  dose = doses,
  b1 = 1,
  b2 = 4,
  b3 = log(5),
  b4 = 5,
  sigma = 0
); d1$sig_emax <- dat$response

dat <- dreamer_data_emax(
  n_cohorts = rep(1, length(doses)),
  dose = doses,
  b1 = 1,
  b2 = 4,
  b3 = log(5),
  b4 = 1,
  sigma = 0
); d1$hyp_emax <- dat$response

dat <- dreamer_data_quad(
  n_cohorts = rep(1, length(doses)),
  dose = doses,
  b1 = 1.3,
  b2 = 25/25,
  b3 = -2/25,
  sigma = 0
); d1$quad <- dat$response

dat <- dreamer_data_logquad(
  n_cohorts = rep(1, length(doses)),
  dose = doses,
  b1 = 1.4,
  b2 = 3*log(5.5)/2,
  b3 = -3/4,
  sigma = 0
); d1$log_quad <- dat$response

dat <- dreamer_data_exp(
  n_cohorts = rep(1, length(doses)),
  dose = doses,
  b1 = 1.1,
  b2 = 2.9,
  b3 = 1,
  sigma = 0
); d1$exp <- dat$response


setEPS(width = 7, height = 5)
devpostscript("figures/DRplot2.eps")
  matplot(d1$dose,d1[,-1],type="l",lwd=2, col=c(1:4,6),
          xlab="dose", ylab="response")
  leg_txt=c("Sigmoid EMAX","Hyperbolic EMAX","Quadratic",
            "Log quadratic","Exponential")
  legend("bottomright",legend=leg_txt,col=c(1:4,6),lty=1:5,lwd=2,
         cex=1,bg="transparent")
dev.off()

# Figure 2
# function 1 ---> ITP
f1 <- function(x, xmax, b1){
  (1 - exp(-b1 * x))/(1 - exp(-b1 * xmax))
}

tmax <- 52 #maximum time T
tm <- seq(0,tmax,0.01)

beta1 <- 0.1
f1_val1 <- sapply(tm,FUN = f1,xmax = tmax, b1 = beta1)
plot(tm, f1_val1, type = "l", lwd=2, col=1, lty=1,
     xlab="time",ylab=expression(f(t)))
beta1 <- 0.25
f1_val2 <- sapply(tm,FUN = f1,xmax = tmax, b1 = beta1)
points(tm, f1_val2, type = "l", lwd=2, col=2, lty=2)
beta1 <- 0.5
f1_val3 <- sapply(tm,FUN = f1,xmax = tmax, b1 = beta1)
points(tm, f1_val3, type = "l", lwd=2, col=3, lty=3)
leg_txt <- c(expression(paste(beta," = 0.10")),
             expression(paste(beta," = 0.25")),
             expression(paste(beta," = 0.50")))
legend("bottomright",legend=leg_txt,col=1:3,lty=1:3,lwd=2)

# Figure 3
# function 1 ---> IDP
fun3.1 <- function(x,d1,d2,b1,b2,gam){
  ((1 - exp(-b1 * x))/(1 - exp(-b1 * d1)))*(x<d1) +
    (1 - gam * ((1 - exp(-b2 * (x - d1)))/(1 - exp(-b2 * (d2 - d1)))))*((x>=d1) & (x<=d2)) +
    (1 - gam)*(x>d2)
}

d1 <- 15
d2 <- 35
b1 <- 0.4
b2 <- -0.1
gam <- 0.05
y3.1 <- fun3.1(tm,b1 = b1,b2=b2, d1=d1, d2=d2, gam = gam)
plot(tm,y3.1,type="l",lwd=2,xlab="time",ylab="function value")
d1 <- 15
d2 <- 35
b1 <- 0.5
b2 <- -0.06
gam <- 0.10
y3.1 <- fun3.1(tm,b1 = b1,b2=b2, d1=d1, d2=d2, gam = gam)
points(tm,y3.1,type="l",lwd=2,col=2,lty=2)
d1 <- 20
d2 <- 40
b1 <- 0.3
b2 <- -0.2
gam <- 0.15
y3.1 <- fun3.1(tm,b1 = b1,b2=b2, d1=d1, d2=d2, gam = gam)
points(tm,y3.1,type="l",lwd=2,col=3,lty=3)
d1 <- 25
d2 <- 45
b1 <- 0.6
b2 <- -0.09
gam <- 0.20
y3.1 <- fun3.1(tm,b1 = b1,b2=b2, d1=d1, d2=d2, gam = gam)
points(tm,y3.1,type="l",lwd=2,col=4,lty=4)
lg <- c(expression(paste(T[1]," = 15; ",T[2]," = 35; ",beta[1]," = 0.4; ",beta[2]," = -0.10; ",gamma," = 0.05")),
        expression(paste(T[1]," = 15; ",T[2]," = 35; ",beta[1]," = 0.5; ",beta[2]," = -0.06; ",gamma," = 0.10")),
        expression(paste(T[1]," = 20; ",T[2]," = 40; ",beta[1]," = 0.3; ",beta[2]," = -0.20; ",gamma," = 0.15")),
        expression(paste(T[1]," = 25; ",T[2]," = 45; ",beta[1]," = 0.6; ",beta[2]," = -0.09; ",gamma," = 0.20"))
)
legend("bottomright",legend=lg,col=1:4,lty=1:4,lwd=2)

# Figure 6
tar_read(scenarios_plot)

# Figure 7
tar_read(mse_plot)

# Figure 8
tar_read(mse_width_plot)

# Figure 9
tar_read(weights_plot)

# Figure 10
tar_read(ocs_plot)

# Figure 11
tar_read(scenarios_plot_binary)

# Figure 12
tar_read(ocs_plot_binary)
