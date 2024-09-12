library(ggplot2)
library(ggpubr)
library(RColorBrewer)

###########################################################################################################################################################################
# Figure 4

calculate_lrt_as_function_of_phi <- function(t, tau, x, N=1000, h=NULL){
  
  # x : vector of data
  # tau : changepoint of interest
  # t : time point to calculate LRT at
  # N : number of phi values in [0, 1] to calculate LRT
  # h : window size  
  
  n <- length(x)
  
  S0 <- mean(x^2)
  S1_tau <- mean(x[1:tau]^2)
  S1_t <- mean(x[1:t]^2)
  
  if ( t <= tau ){
    a1 <- n * S0 * S1_t / (tau * S1_tau)
    a0 <- 0
  } else {
    a1 <- n * S0 * (n * S0 - t * S1_t) / (t * n * S0 - t * tau * S1_tau)
    a0 <- n * S0 * (t * S1_t - tau * S1_tau) / (t * n * S0 - t * tau * S1_tau)
  }
  
  phi <- seq(1/N, 1 - 1/N, length.out=(N - 1))
  lrt <- rep(NA, N - 1)
  for ( i in 1:length(phi) ){
    lrt[i] <- n * log(S0) - t * log(a1 * phi[i] + a0) - (n - t) * log(n/(n - t) * S0 - t/(n - t) * (a1 * phi[i] + a0))
  }
  
  return(lrt)
}

seeds <- c(1, 4)
n <- 20
N_phi <- 1000
phi <- seq(1/N_phi, 1 - 1/N_phi, length.out=(N_phi - 1))
plots_list <- as.list(rep(NA, 2))

for ( plot_index in 1:2 ){
  
  set.seed(seeds[plot_index])
  
  x <- c(rnorm(n/2), rnorm(n/2, sd=2))
  lrs <- n * log(mean(x^2)) - (1:(n - 1)) * log(cumsum(x[-n]^2)/(1:(n - 1))) - (n - 1):1 * log(cumsum(x[n:2]^2)[(n - 1):1]/((n - 1):1))
  tau <- which.max(lrs)
  phi_obs <- sum(x[1:tau]^2) / sum(x^2)
  
  lrs_grid <- apply(matrix(1:(n - 1), ncol=1), 1, calculate_lrt_as_function_of_phi, x=x, tau=tau, N=N_phi)
  
  dat <- data.frame(phi=rep(phi, n - 1), lrs=c(lrs_grid), tau=as.factor(rep(1:(n - 1), each=N_phi - 1)))
  tau_which_max_lrs <- apply(lrs_grid, 1, which.max)
  possible_taus <- unique(tau_which_max_lrs)
  colours <- brewer.pal(length(possible_taus), "Set1")
  
  plots_list[[plot_index]] <- ggplot(dat)
  for ( k in 1:length(possible_taus) ){
    phi_inds_where_tau_is_cp <- which(tau_which_max_lrs == possible_taus[k])
    num_phi_inds <- length(phi_inds_where_tau_is_cp)
    if ( sum(phi_inds_where_tau_is_cp[-1] - phi_inds_where_tau_is_cp[-num_phi_inds] != 1) == 0 ){
      plots_list[[plot_index]] <- plots_list[[plot_index]] + annotate("rect", xmin=phi[phi_inds_where_tau_is_cp[1]] - 1/(2 * N_phi), 
                                                                      xmax=phi[phi_inds_where_tau_is_cp[num_phi_inds]] + 1/(2 * N_phi), ymin=-Inf, ymax=Inf, fill=colours[k], alpha=0.15)
    } else {
      zz <- c(0, which(phi_inds_where_tau_is_cp[-1] - phi_inds_where_tau_is_cp[-num_phi_inds] != 1), num_phi_inds)
      for ( i in 1:(length(zz) - 1) ){
        plots_list[[plot_index]] <- plots_list[[plot_index]] + annotate("rect", xmin=phi[phi_inds_where_tau_is_cp[zz[i] + 1]] - 1/(2 * N_phi), 
                                                                        xmax=phi[phi_inds_where_tau_is_cp[zz[i + 1]]] + 1/(2 * N_phi), ymin=-Inf, ymax=Inf, fill=colours[k], alpha=0.15)
      }
    }
  }
  plots_list[[plot_index]] <- plots_list[[plot_index]] + geom_line(aes(x=phi, y=lrs, group=tau), colour="grey")
  for ( k in length(possible_taus):1 ){
    plots_list[[plot_index]] <- plots_list[[plot_index]] + geom_line(data=dat[dat$tau == possible_taus[k],], aes(x=phi, y=lrs), colour=colours[k], linewidth=1.1)
  }
  plots_list[[plot_index]] <- plots_list[[plot_index]] + theme_classic() + labs(x="phi", y="") + theme(axis.text=element_text(size=20), axis.title=element_text(size=22))
}

# Plot
ggarrange(plots_list[[1]] + labs(y="Likelihood ratio statistic"), plots_list[[2]])

###########################################################################################################################################################################
# Figure 5: Gaussian process plots

# Generate X
set.seed(29)
x <- c(rnorm(100), rnorm(100, sd=2), rnorm(100))
results <- find_changepoints(x, method="bs", params=list(loss="lrs", threshold=10), model="var")
b <- results$changepoints[1]
h <- 50
N <- 100

h <- 50
phi_obs <- sum(x[(b - h + 1):b]^2) / sum(x[(b - h + 1):(b + h)]^2)
phi_star <- qbeta(1 - pbeta(phi_obs, h/2, h/2), h/2, h/2)
phi_upper <- max(phi_obs, phi_star)
phi_lower <- min(phi_obs, phi_star)

# Generate phi and calculate I(phi in S)
N <- 100
phi <- runif(N) / N + seq(0, 1 - 1/N, length.out=N)
alpha <- sqrt(phi/phi_obs)
beta <- sqrt((1 - phi)/(1 - phi_obs))
inS <- rep(NA, N)
for ( i in 1:N ){
  x2 <- x
  x2[(b - h + 1):(b + h)] <- c(alpha[i] * x[(b - h + 1):b], beta[i] * x[(b + 1):(b + h)])
  results_phi <- find_changepoints(x2, method="bs", params=list(threshold=results$params$threshold, maxiter=results$params$maxiter, loss="lrs"), model="var")
  if ( length(results_phi$changepoints) >= 1 ){
    inS[i] <- ifelse(b %in% results_phi$changepoints, 1, 0)
  } else {
    inS[i] <- 0
  }
}

g1 <- ggplot(data.frame(phi=phi, inS=inS)) + 
  geom_point(aes(x=phi, y=inS), colour="red") + 
  geom_line(aes(x=phi, y=inS)) +
  theme_classic() + labs(x="phi", y="I(phi in S)") + 
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

# Implement GP
phi2 <- sort(c(phi, runif(10000)))
posterior_results <- calculate_posterior(phi2, phi, inS, l=1)
mu_pos <- posterior_results$mu
cov_pos <- posterior_results$cov

g2 <- ggplot(data.frame(phi2=phi2, mu=mu_pos)) +  
  geom_line(aes(x=phi2, y=mu + qnorm(0.975) * diag(cov_pos)), colour="red") +
  geom_line(aes(x=phi2, y=mu + qnorm(0.025) * diag(cov_pos)), colour="red") +
  geom_line(aes(x=phi2, y=mu), colour="blue") +
  theme_classic() + labs(x="phi", y="p_hat") + 
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

# Estimate q(phi)
q2 <- mu_pos * dbeta(phi2, h/2, h/2)
q2[q2 < 0] <- 0
q2 <- q2 / sum(q2) * length(q2)

g3 <- ggplot(data.frame(phi2, q2, pi=dbeta(phi2, h/2, h/2))) + 
  geom_line(aes(x=phi2, y=q2)) + 
  geom_line(aes(x=phi2, y=pi), colour="red", linetype="dotted") +
  theme_classic() + labs(x="phi", y="q_hat") + 
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

ggarrange(g1, g2, g3, ncol=3)
