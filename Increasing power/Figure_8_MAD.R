# QQ plots where the variance is estimated using MAD

###################################################################################################################
# Fig. 7(a), H0 with h = 10

set.seed(100)

h <- 10
n <- 1000
num_iter <- 1000
pvals <- matrix(NA, nrow=num_iter, ncol=6)

iter <- 1
while( iter <= num_iter ){
  x <- rnorm(n)
  results <- find_changepoints(x, "bs", list(threshold=3, maxiter=1))
  if ( length(results$changepoints) >= 1 ){
    if ( results$changepoints[1] >= h & results$changepoints[1] <= (n - h) ){
      sig2 <- mad(x)
      results_pvals <- calculate_pvals_all(results, N=50, h=h, sigma2=sig2, num_pvals=1, include_original=TRUE, return_probs=TRUE)
      pvals[iter, 1] <- results_pvals$P_both[1] / results_pvals$P_phi_in_S[1]
      pvals[iter, 2] <- sum(results_pvals$P_both[1:2]) / sum(results_pvals$P_phi_in_S[1:2])
      pvals[iter, 3] <- sum(results_pvals$P_both[1:5]) / sum(results_pvals$P_phi_in_S[1:5])
      pvals[iter, 4] <- sum(results_pvals$P_both[1:10]) / sum(results_pvals$P_phi_in_S[1:10])
      pvals[iter, 5] <- sum(results_pvals$P_both[1:20]) / sum(results_pvals$P_phi_in_S[1:20])
      pvals[iter, 6] <- sum(results_pvals$P_both) / sum(results_pvals$P_phi_in_S)
      iter <- iter + 1
      print(iter)
    }
  }
}
pvals[pvals > 1] <- 1 # some p-values are slightly above 1 because of floating point errors
pvals <- pvals[rowSums(is.na(pvals)) == 0, ]
colnames(pvals) <- paste0("N", c(1, 2, 5, 10, 20, 50))
saveRDS(pvals, paste0("H0_inc_orig_bs_n", n, "_h", h, "_N", num_iter, "_mad.rds"))

###################################################################################################################
# under H1

set.seed(100)

n <- 1000

for ( delta in c(1, 2, 3) ){
  for ( h in c(10, 20, 50) ){
    pvals <- matrix(NA, nrow=1000, ncol=6)
    iter <- 1
    while( iter <= n ){
      x <- rnorm(n) + c(rep(0, n/2), rep(delta, n/2))
      results <- find_changepoints(x, "bs", list(threshold=3, maxiter=1))
      if ( length(results$changepoints) >= 1 ){
        if ( results$changepoints[1] >= h & results$changepoints[1] <= n - h ){
          sigma2 <- mad(x)
          pval_results <- calculate_pvals_all(results, h=h, sigma2=sigma2, N=50, include_original=TRUE, num_pvals=1, return_probs=TRUE)
          pvals[iter, 1] <- pval_results$P_both[1] / pval_results$P_phi_in_S[1]
          pvals[iter, 2] <- sum(pval_results$P_both[1:2]) / sum(pval_results$P_phi_in_S[1:2])
          pvals[iter, 3] <- sum(pval_results$P_both[1:5]) / sum(pval_results$P_phi_in_S[1:5])
          pvals[iter, 4] <- sum(pval_results$P_both[1:10]) / sum(pval_results$P_phi_in_S[1:10])
          pvals[iter, 5] <- sum(pval_results$P_both[1:20]) / sum(pval_results$P_phi_in_S[1:20])
          pvals[iter, 6] <- sum(pval_results$P_both) / sum(pval_results$P_phi_in_S)
          iter <- iter + 1
          print(iter)
      }
    }
    pvals[pvals > 1] <- 1 # some p-values are slightly above 1 because of floating point errors
    colnames(pvals) <- paste0("N", c(1, 2, 5, 10, 20, 50))
    saveRDS(pvals, paste0("pvals_T", n, "_bs_delta", delta, "_h", h, "_mad.rds"))
  }
}


###################################################################################################################
# Fig. 7(a), H0 with h = 10

pvals <- readRDS(paste0(pvals_dir, "/H0_inc_orig_bs_n1000_h10_N1000_mad.rds"))

# Sort p-values
pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=5)
for ( i in 1:5 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

# Plot p-values against theoretical quantities from U(0,1)
dat <- data.frame(z=c(rep(seq(0, 1, length.out=nrow(pvals_sorted)), 5), NA), p=c(as.vector(pvals_sorted), NA), 
                  N=as.factor(c(rep(c(1, 5, 10, 20, 50), each=nrow(pvals_sorted)), 2)))

fig7a <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="P-value") + 
  theme_classic() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.title=element_text(size=16), legend.text=element_text(size=16))

###################################################################################################################
# Fig. 7(b), H1 with delta = 2 and h = 10

h <- 10
delta <- 2
method <- "bs"
pvals <- readRDS(paste0(pvals_dir, "/pvals_T1000_bs_delta", delta, "_h", h, "_mad.rds"))
pvals_sorted <- matrix(NA, nrow=1000, ncol=6)
for ( i in 1:6 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}
dat <- data.frame(p1=c(rep(pvals_sorted[,1], 5), NA), p2=c(as.vector(pvals_sorted[,2:6]), NA), 
                  N=as.factor(c(rep(c(2, 5, 10, 20, 50), each=1000), 1)))

fig7b <- ggplot(dat) + geom_line(aes(x=p1, y=p2, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="P-value (N = 1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18))

###################################################################################################################
# Combine plots

pdf(paste0(plots_dir, "/fig7_qq_estimated_variance.pdf"), width=10, height=6)
ggarrange(fig7a + theme(legend.key.size=unit(3, "line")) , fig7b, ncol=2, common.legend=TRUE, legend="bottom")
dev.off()