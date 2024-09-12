library(ggplot2)
library(ggpubr)

kernel <- function(X1, X2, l=1, sigma_f=1){
  dists <- abs(outer(X1, X2, "-"))
  return(sigma_f^2 * exp(-0.5 / l^2 * dists))
}

posterior <- function(X_s, X_train, Y_train, l=1, sigma_f=1, sigma_y=10^(-8)){
  
  K <- kernel(X_train, X_train, l, sigma_f) + sigma_y ^ 2 * diag(length(X_train))
  K_s <- kernel(X_train, X_s, l, sigma_f)
  K_ss <- kernel(X_s, X_s, l, sigma_f) + 10^(-8) * diag(length(X_s))
  K_inv <- solve(K)
  
  mu_s <- t(K_s) %*% K_inv %*% Y_train
  cov_s <- K_ss - t(K_s) %*% K_inv %*% K_s
  
  return(list(mu_s=mu_s, cov_s=cov_s))
}

N <- 100
h <- 50

set.seed(29)
x <- c(rnorm(100), rnorm(100, sd=2), rnorm(100))
results <- find_changepoints(x, "bs", list(threshold=10, loss="lrs"), model="var")
b <- results$b[1]
    
phi_obs <- sum(x[(b - h + 1):b]^2) / sum(x[(b - h + 1):(b + h)]^2)
phi_star <- qbeta(1 - pbeta(phi_obs, h/2, h/2), h/2, h/2)
phi_upper <- max(phi_obs, phi_star)
phi_lower <- min(phi_obs, phi_star)

phi <- runif(N) / N + seq(0, (N - 1)/N, length.out=N)
alpha <- sqrt(phi/phi_obs)
beta <- sqrt((1 - phi)/(1 - phi_obs))
inS <- rep(NA, N)
for ( i in 1:N ){
  x2 <- x
  x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
  r2 <- find_cps_var(x2, method="bs", c(results$params$threshold, results$params$maxiter), loss="lrs")
  if ( length(r2$cps) >= 1 ){
    inS[i] <- ifelse(b %in% r2$cps, 1, 0)
  } else {
    inS[i] <- 0
  }
}
    
X <- sort(c(phi, runif(1000)))
l_list <- c(0.01, 0.1, 1, 10, 100, 1000)
g1_list <- g2_list <- as.list(rep(0, length(l_list)))
for ( i in 1:length(l_list) ){
  posterior_results <- posterior(X, phi, inS, l=l_list[i])
  mu_star <- posterior_results$mu_s
  mu_star[mu_star < 0] <- 0
  mu_star[mu_star > 1] <- 1
  pi_hat <- mu_star * dbeta(X, h/2, h/2)
  g1_list[[i]] <- ggplot(data.frame(phi=X, mu=mu_star)) + geom_line(aes(x=phi, y=mu)) + 
                    geom_point(data=data.frame(phi=phi, y=inS), aes(x=phi, y=y), colour="red") +
                    theme_classic() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
}


ggarrange(g1_list[[1]] + labs(x="", y="p_hat"), g1_list[[2]] + labs(x="", y=""), g1_list[[3]] + labs(x="", y=""),
          g1_list[[4]] + labs(x="phi", y="p_hat"), g1_list[[5]] + labs(x="phi", y=""), g1_list[[6]] + labs(x="phi", y=""))
