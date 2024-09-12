# Figure 2: QQ plots for binary segmentation under H0

###################################################################################################################

# Simulating data under H0, estimating a change, and calculating p-values with and without the observed value of psi being included

set.seed(100)

n <- 1000 # size of dataset
num_iter <- 1000 # number of simulations

pvals_not_10 <- pvals_not_20 <- pvals_obs_10 <- pvals_obs_20 <- matrix(NA, nrow=num_iter, ncol=5)
colnames(pvals_not_10) <- colnames(pvals_not_20) <- colnames(pvals_obs_10) <- colnames(pvals_obs_20) <- paste0("N", c(1, 5, 10, 20, 50))

h <- 10 # window size
iter <- 1
while( iter <= num_iter ){
  x <- rnorm(n)
  results <- find_changepoints(x, "bs", list(threshold=2, maxiter=1))
  if ( length(results$changepoints) >= 1 ){ # only keep simulations for which the changepoint algorithm detects at least one changepoint
    if ( results$changepoints[1] >= h & results$changepoints[1] <= n - h ){
      pval_results <- calculate_pvals_all(results, N=51, h=h, sigma2=1, num_pvals=1, include_original=TRUE, return_probs=TRUE)
      
      # Not including observed phi
      pvals_not_10[iter, 1] <- pval_results$P_both[2] / pval_results$P_phi_in_S[2]
      pvals_not_10[iter, 2] <- sum(pval_results$P_both[2:6]) / sum(pval_results$P_phi_in_S[2:6])
      pvals_not_10[iter, 3] <- sum(pval_results$P_both[2:11]) / sum(pval_results$P_phi_in_S[2:11])
      pvals_not_10[iter, 4] <- sum(pval_results$P_both[2:21]) / sum(pval_results$P_phi_in_S[2:21])
      pvals_not_10[iter, 5] <- sum(pval_results$P_both[2:51]) / sum(pval_results$P_phi_in_S[2:51])
      
      # Including observed phi
      pvals_obs_10[iter, 1] <- pval_results$P_both[1] / pval_results$P_phi_in_S[1]
      pvals_obs_10[iter, 2] <- sum(pval_results$P_both[1:5]) / sum(pval_results$P_phi_in_S[1:5])
      pvals_obs_10[iter, 3] <- sum(pval_results$P_both[1:10]) / sum(pval_results$P_phi_in_S[1:10])
      pvals_obs_10[iter, 4] <- sum(pval_results$P_both[1:20]) / sum(pval_results$P_phi_in_S[1:20])
      pvals_obs_10[iter, 5] <- sum(pval_results$P_both[1:50]) / sum(pval_results$P_phi_in_S[1:50])
      
      print(iter)
      iter <- iter + 1
    }
  }
}
pvals_not_10[pvals_not_10 > 1] <- 1 # some p-values may be slightly above 1 due to floating point errors
pvals_obs_10[pvals_obs_10 > 1] <- 1



h <- 20 # window size
iter <- 1
while( iter <= num_iter ){
  x <- rnorm(n)
  results <- find_changepoints(x, "bs", list(threshold=2, maxiter=1))
  if ( length(results$changepoints) >= 1 ){ # only keep simulations for which the changepoint algorithm detects at least one changepoint
    if ( results$changepoints[1] >= h & results$changepoints[1] <= n - h ){
      pval_results <- calculate_pvals_all(results, N=51, h=h, sigma2=1, num_pvals=1, include_original=TRUE, return_probs=TRUE)
      
      # Not including observed phi
      pvals_not_20[iter, 1] <- pval_results$P_both[2] / pval_results$P_phi_in_S[2]
      pvals_not_20[iter, 2] <- sum(pval_results$P_both[2:6]) / sum(pval_results$P_phi_in_S[2:6])
      pvals_not_20[iter, 3] <- sum(pval_results$P_both[2:11]) / sum(pval_results$P_phi_in_S[2:11])
      pvals_not_20[iter, 4] <- sum(pval_results$P_both[2:21]) / sum(pval_results$P_phi_in_S[2:21])
      pvals_not_20[iter, 5] <- sum(pval_results$P_both[2:51]) / sum(pval_results$P_phi_in_S[2:51])
      
      # Including observed phi
      pvals_obs_20[iter, 1] <- pval_results$P_both[1] / pval_results$P_phi_in_S[1]
      pvals_obs_20[iter, 2] <- sum(pval_results$P_both[1:5]) / sum(pval_results$P_phi_in_S[1:5])
      pvals_obs_20[iter, 3] <- sum(pval_results$P_both[1:10]) / sum(pval_results$P_phi_in_S[1:10])
      pvals_obs_20[iter, 4] <- sum(pval_results$P_both[1:20]) / sum(pval_results$P_phi_in_S[1:20])
      pvals_obs_20[iter, 5] <- sum(pval_results$P_both[1:50]) / sum(pval_results$P_phi_in_S[1:50])
      
      print(iter)
      iter <- iter + 1
    }
  }
}
pvals_not_20[pvals_not_20 > 1] <- 1 # some p-values may be slightly above 1 due to floating point errors
pvals_obs_20[pvals_obs_20 > 1] <- 1

###################################################################################################################
# Figure 2(a): h = 10, not including original

# Create dataframe with sorted p-values
num_pvals <- num_iter - colSums(is.na(pvals_not_10))
dat <- data.frame(p=sort(pvals_not_10[,1]), z=seq(0, 1, length.out=num_pvals[1]))
for ( i in 2:5 ){
  dat <- rbind(dat, data.frame(p=sort(pvals_not_10[,i]), z=seq(0, 1, length.out=num_pvals[i])))
}
dat$N <- as.factor(rep(c(1, 5, 10, 20, 50), num_pvals))

# Plot p-values against theoretical quantities from U(0,1)
fig2a <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="P-value") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

###################################################################################################################
# Figure 2(b): h = 20, not including original

# Create dataframe with sorted p-values
num_pvals <- num_iter - colSums(is.na(pvals_not_20))
dat <- data.frame(p=sort(pvals_not_20[,1]), z=seq(0, 1, length.out=num_pvals[1]))
for ( i in 2:5 ){
  dat <- rbind(dat, data.frame(p=sort(pvals_not_20[,i]), z=seq(0, 1, length.out=num_pvals[i])))
}
dat$N <- as.factor(rep(c(1, 5, 10, 20, 50), num_pvals))

# Plot p-values against theoretical quantities from U(0,1)
fig2b <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="P-value") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20), legend.title=element_text(size=16), legend.text=element_text(size=16))

###################################################################################################################
# Figure 2(c): h = 10, including observed

# Create dataframe with sorted p-values
num_pvals <- num_iter - colSums(is.na(pvals_obs_10))
dat <- data.frame(p=sort(pvals_obs_10[,1]), z=seq(0, 1, length.out=num_pvals[1]))
for ( i in 2:5 ){
  dat <- rbind(dat, data.frame(p=sort(pvals_obs_10[,i]), z=seq(0, 1, length.out=num_pvals[i])))
}
dat$N <- as.factor(rep(c(1, 5, 10, 20, 50), num_pvals))

# Plot p-values against theoretical quantities from U(0,1)
fig2c <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

###################################################################################################################
# Figure 2(d): h = 20, including observed

# Create dataframe with sorted p-values
num_pvals <- num_iter - colSums(is.na(pvals_obs_20))
dat <- data.frame(p=sort(pvals_obs_20[,1]), z=seq(0, 1, length.out=num_pvals[1]))
for ( i in 2:5 ){
  dat <- rbind(dat, data.frame(p=sort(pvals_obs_20[,i]), z=seq(0, 1, length.out=num_pvals[i])))
}
dat$N <- as.factor(rep(c(1, 5, 10, 20, 50), num_pvals))

# Plot p-values against theoretical quantities from U(0,1)
fig2d <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))


###################################################################################################################
# Plot Figure 2

ggarrange(fig2a + theme(legend.title=element_text(size=16), legend.text=element_text(size=16), legend.key.size=unit(3, "line")), 
          fig2b, fig2c, fig2d, ncol=4, nrow=1, common.legend=TRUE, legend="bottom")