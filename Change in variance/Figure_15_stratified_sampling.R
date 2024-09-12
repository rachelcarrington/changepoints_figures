library(ggplot2)
library(ggpubr)

set.seed(200)
h <- 20
threshold <- 10
maxiter <- 1
N <- 200
NN <- 1000

x <- c(rnorm(100), rnorm(100, sd=2))
results <- find_changepoints(x, "bs", list(threshold=threshold, maxiter=maxiter, loss="lrs"), model="var")
b <- results$b[1]
phi_obs <- sum(x[(b - h + 1):b]^2) / sum(x[(b - h + 1):(b + h)]^2)
phi_star <- qbeta(1 - pbeta(phi_obs, h/2, h/2), h/2, h/2)
phi_upper <- max(phi_obs, phi_star)
phi_lower <- min(phi_obs, phi_star)

p_hat <- matrix(NA, nrow=NN, ncol=3)

for ( iter in 1:NN ){
  phi1 <- rbeta(N, h/2, h/2)
  alpha <- sqrt(phi1/phi_obs)
  beta <- sqrt((1 - phi1)/(1 - phi_obs))
  inS1 <- rep(NA, N)
  for ( i in 1:N ){
    x2 <- x
    x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
    r2 <- find_cps_var(x2, method="bs", c(results$params$threshold, results$params$maxiter), loss="lrs")
    if ( length(r2$cps) >= 1 ){
      inS1[i] <- ifelse(b %in% r2$cps, 1, 0)
    } else {
      inS1[i] <- 0
    }
  }
  p_hat[iter, 1] <- sum(inS1 == 1 & !(phi1 > phi_lower & phi1 < phi_upper)) / sum(inS1 == 1)

  phi2 <- rbeta(N, h/10, h/10)
  alpha <- sqrt(phi2/phi_obs)
  beta <- sqrt((1 - phi2)/(1 - phi_obs))
  inS2 <- rep(NA, N)
  for ( i in 1:N ){
    x2 <- x
    x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
    r2 <- find_cps_var(x2, method="bs", c(results$params$threshold, results$params$maxiter), loss="lrs")
    if ( length(r2$cps) >= 1 ){
      inS2[i] <- ifelse(b %in% r2$cps, 1, 0)
    } else {
      inS2[i] <- 0
    }
  }
  top <- sum(dbeta(phi2[inS2 == 1 & !(phi2 > phi_lower & phi2 < phi_upper)], h/2, h/2) / dbeta(phi2[inS2 == 1 & !(phi2 > phi_lower & phi2 < phi_upper)], h/10, h/10))
  bottom <- sum(dbeta(phi2[inS2 == 1], h/2, h/2) / dbeta(phi2[inS2 == 1], h/10, h/10))
  p_hat[iter, 2] <- top / bottom

  phi3 <- rbeta(N, h/100, h/100)
  alpha <- sqrt(phi3/phi_obs)
  beta <- sqrt((1 - phi3)/(1 - phi_obs))
  inS3 <- rep(NA, N)
  for ( i in 1:N ){
    x2 <- x
    x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
    r2 <- find_cps_var(x2, method="bs", c(results$params$threshold, results$params$maxiter), loss="lrs")
    if ( length(r2$cps) >= 1 ){
      inS3[i] <- ifelse(b %in% r2$cps, 1, 0)
    } else {
      inS3[i] <- 0
    }
  }
  top <- sum(dbeta(phi3[inS3 == 1 & !(phi3 > phi_lower & phi3 < phi_upper)], h/2, h/2) / dbeta(phi3[inS3 == 1 & !(phi3 > phi_lower & phi3 < phi_upper)], h/100, h/100))
  bottom <- sum(dbeta(phi3[inS3 == 1], h/2, h/2) / dbeta(phi3[inS3 == 1], h/100, h/100))
  p_hat[iter, 3] <- top / bottom
  print(iter)
}

dat <- data.frame(z=c(seq(0, N, length.out=sum(!is.na(p_hat[,1]))), seq(0, N, length.out=sum(!is.na(p_hat[,2]))), seq(0, N, length.out=sum(!is.na(p_hat[,3])))),
                  p=c(sort(p_hat[,1]), sort(p_hat[,2]), sort(p_hat[,3])),
                  k=as.factor(c(rep(1, sum(!is.na(p_hat[,1]))), rep(5, sum(!is.na(p_hat[,2]))), rep(50, sum(!is.na(p_hat[,3]))))))

g200 <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=k, linetype=k), linewidth=1.2) + 
    geom_hline(yintercept=mean(p_hat, na.rm=TRUE), colour="grey", linetype="dashed") +
    theme_classic() + labs(x="Index", y="P-value estimate") +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=20), legend.text=element_text(size=18), legend.title=element_text(size=18),
          legend.position="bottom", legend.key.size=unit(3, "line"))


## stratified sampling

p_hat <- matrix(NA, nrow=NN, ncol=3)

for ( iter in 1:NN ){
  z1 <- runif(N)/N + seq(0, 1 - 1/N, by=1/N)
  phi1 <- qbeta(z1, h/2, h/2)
  alpha <- sqrt(phi1/phi_obs)
  beta <- sqrt((1 - phi1)/(1 - phi_obs))
  inS1 <- rep(NA, N)
  for ( i in 1:N ){
    x2 <- x
    x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
    r2 <- find_cps_var(x2, method="bs", c(results$params$threshold, results$params$maxiter), loss="lrs")
    if ( length(r2$cps) >= 1 ){
      inS1[i] <- ifelse(b %in% r2$cps, 1, 0)
    } else {
      inS1[i] <- 0
    }
  }
  p_hat[iter, 1] <- sum(inS1 == 1 & !(phi1 > phi_lower & phi1 < phi_upper)) / sum(inS1 == 1)

  z2 <- runif(N)/N + seq(0, 1 - 1/N, by=1/N)
  phi2 <- qbeta(z2, h/10, h/10)
  alpha <- sqrt(phi2/phi_obs)
  beta <- sqrt((1 - phi2)/(1 - phi_obs))
  inS2 <- rep(NA, N)
  for ( i in 1:N ){
    x2 <- x
    x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
    r2 <- find_cps_var(x2, method="bs", c(results$params$threshold, results$params$maxiter), loss="lrs")
    if ( length(r2$cps) >= 1 ){
      inS2[i] <- ifelse(b %in% r2$cps, 1, 0)
    } else {
      inS2[i] <- 0
    }
  }
  top <- sum(dbeta(phi2[inS2 == 1 & !(phi2 > phi_lower & phi2 < phi_upper)], h/2, h/2) / dbeta(phi2[inS2 == 1 & !(phi2 > phi_lower & phi2 < phi_upper)], h/10, h/10))
  bottom <- sum(dbeta(phi2[inS2 == 1], h/2, h/2) / dbeta(phi2[inS2 == 1], h/10, h/10))
  p_hat[iter, 2] <- top / bottom

  z3 <- runif(N)/N + seq(0, 1 - 1/N, by=1/N)
  phi3 <- qbeta(z3, h/100, h/100)
  alpha <- sqrt(phi3/phi_obs)
  beta <- sqrt((1 - phi3)/(1 - phi_obs))
  inS3 <- rep(NA, N)
  for ( i in 1:N ){
    x2 <- x
    x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
    r2 <- find_cps_var(x2, method="bs", c(results$params$threshold, results$params$maxiter), loss="lrs")
    if ( length(r2$cps) >= 1 ){
      inS3[i] <- ifelse(b %in% r2$cps, 1, 0)
    } else {
      inS3[i] <- 0
    }
  }
  top <- sum(dbeta(phi3[inS3 == 1 & !(phi3 > phi_lower & phi3 < phi_upper)], h/2, h/2) / dbeta(phi3[inS3 == 1 & !(phi3 > phi_lower & phi3 < phi_upper)], h/100, h/100))
  bottom <- sum(dbeta(phi3[inS3 == 1], h/2, h/2) / dbeta(phi3[inS3 == 1], h/100, h/100))
  p_hat[iter, 3] <- top / bottom
  print(iter)
}


dat <- data.frame(z=c(seq(0, N, length.out=sum(!is.na(p_hat[,1]))), seq(0, N, length.out=sum(!is.na(p_hat[,2]))), seq(0, N, length.out=sum(!is.na(p_hat[,3])))),
                  p=c(sort(p_hat[,1]), sort(p_hat[,2]), sort(p_hat[,3])),
                  k=as.factor(c(rep(1, sum(!is.na(p_hat[,1]))), rep(5, sum(!is.na(p_hat[,2]))), rep(50, sum(!is.na(p_hat[,3]))))))

g_strat <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=k, linetype=k), linewidth=1.2) + 
    geom_hline(yintercept=mean(p_hat, na.rm=TRUE), colour="grey", linetype="dashed") +
    theme_classic() + labs(x="Index", y="") +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=20), legend.text=element_text(size=18), legend.title=element_text(size=18),
          legend.position="bottom", legend.key.size=unit(3, "line"))

ggarrange(g200 + coord_cartesian(ylim=c(0, 1)), g_strat + coord_cartesian(ylim=c(0, 1)), ncol=2, common.legend=TRUE, legend="bottom")
