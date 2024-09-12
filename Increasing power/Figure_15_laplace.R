# Figure 15: QQ plots, Laplace distribution

###################################################################################################################
# H0 

set.seed(100)

h_list <- c(10, 20)
s_list <- c(0.5, 1, 2)
n <- 200 # size of dataset
num_iter <- 1000 # number of simulations

for ( h in h_list ){
  for ( s in s_list ){
    pvals <- matrix(NA, nrow=num_iter, ncol=5)
    iter <- 1
    while( iter <= num_iter ){
      x <- rlaplace(n, s=s)
      results <- find_changepoints(x, "bs", list(threshold=4 * s, maxiter=1))
      if ( length(results$changepoints) >= 1 ){ # only keep simulations for which the changepoint algorithm detects at least one changepoint
        if ( results$changepoints[1] >= h & results$changepoints[1] <= (n - h) ){
          pval_results <- calculate_pvals_all(results, N=50, h=h, sigma2=2*s^2, num_pvals=1, include_original=TRUE, return_probs=TRUE)
          pvals[iter,1] <- pval_results$P_both[1] / pval_results$P_phi_in_S[1]
          pvals[iter,2] <- sum(pval_results$P_both[1:5]) / sum(pval_results$P_phi_in_S[1:5])
          pvals[iter,3] <- sum(pval_results$P_both[1:10]) / sum(pval_results$P_phi_in_S[1:10])
          pvals[iter,4] <- sum(pval_results$P_both[1:20]) / sum(pval_results$P_phi_in_S[1:20])
          pvals[iter,5] <- sum(pval_results$P_both) / sum(pval_results$P_phi_in_S)
          iter <- iter + 1
          print(iter)
        }
      }
    }
    pvals[pvals > 1] <- 1 # some p-values are slightly above 1 because of floating point errors
    pvals <- pvals[rowSums(is.na(pvals)) == 0,]
    colnames(pvals) <- paste0("N", c(1, 5, 10, 20, 50))
    write.csv(pvals, paste0("pvals_laplace_H0_T", n, "_h", h, "_s", s, ".csv"))
  }
}

###################################################################################################################
# H1

set.seed(100)

for ( delta in c(1, 2) ){
  for ( h in c(10, 20) ){
    for ( s in c(0.5, 1, 2) ){
      pvals <- matrix(NA, nrow=1000, ncol=6)
      iter <- 1
      while( iter <= 1000 ){
        x <- rlaplace(200, s=s) + c(rep(delta/2, 100), rep(-delta/2, 100))
        results <- find_changepoints(x, "bs", list(threshold=2 * s, maxiter=1))
        if ( length(results$changepoints) >=1 ){
          if ( results$changepoints[1] >= h & results$changepoints[1] <= (200 - h) ){
            pval_results <- calculate_pvals_all(results, N=50, h=h, sigma2=2*s^2, num_pvals=1, include_original=TRUE, return_probs=TRUE)
            pvals[iter,1] <- pval_results$P_both[1] / pval_results$P_phi_in_S[1]
            pvals[iter,2] <- sum(pval_results$P_both[1:2]) / sum(pval_results$P_phi_in_S[1:2])
            pvals[iter,3] <- sum(pval_results$P_both[1:5]) / sum(pval_results$P_phi_in_S[1:5])
            pvals[iter,4] <- sum(pval_results$P_both[1:10]) / sum(pval_results$P_phi_in_S[1:10])
            pvals[iter,5] <- sum(pval_results$P_both[1:20]) / sum(pval_results$P_phi_in_S[1:20])
            pvals[iter,6] <- sum(pval_results$P_both) / sum(pval_results$P_phi_in_S)
            iter <- iter + 1
            print(iter)
          }
        }
      }
      pvals[pvals > 1] <- 1 # some p-values are slightly above 1 because of floating point errors
      colnames(pvals) <- paste0("N", c(1, 2, 5, 10, 20, 50))
      
      write.csv(pvals, paste0("pvals_T200_laplace_delta", delta, "_h", h, "_s", s, ".csv")) 
    }
  }
}


###################################################################################################################
# Plots

# Fig. 15(a): H0, h = 10, s = 0.5
h <- 10
s <- 0.5
n <- 200
pvals <- read.csv(paste0(pvals_dir, "/pvals_laplace_H0_T", n, "_h", h, "_s", s, ".csv"))[,-1]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:5 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

# Plot p-values against theoretical quantities from U(0,1)
dat1 <- data.frame(x=rep(seq(0, 1, length.out=nrow(pvals_sorted)), 4), y=as.vector(pvals_sorted[,c(1,2,3,5)]), 
                   N=rep(as.factor(c(1, 5, 10, 50)), each=nrow(pvals_sorted)))

fig14a <- ggplot(dat1) + geom_abline(intercept=0, slope=1) + 
  geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  labs(x="U(0,1)", y="P-value", colour="Number of samples", linetype="Number of samples") + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22), legend.text=element_text(size=22), legend.title=element_text(size=22),
        legend.key.size=unit(4, "line"))

# Fig. 15(b): H0, h = 10, s = 1
h <- 10
s <- 1
n <- 200
pvals <- read.csv(paste0(pvals_dir, "/pvals_laplace_H0_T", n, "_h", h, "_s", s, ".csv"))[,-1]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:5 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

# Plot p-values against theoretical quantities from U(0,1)
dat2 <- data.frame(x=rep(seq(0, 1, length.out=nrow(pvals_sorted)), 4), y=as.vector(pvals_sorted[,c(1,2,3,5)]), 
                   N=rep(as.factor(c(1, 5, 10, 50)), each=nrow(pvals_sorted)))

fig14b <- ggplot(dat2) + geom_abline(intercept=0, slope=1) + 
  geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  labs(x="U(0,1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))


# Fig. 15(c): H1, h = 10, s = 0.5, delta = 1
h <- 10
s <- 0.5
n <- 200
delta <- 1
pvals <- read.csv(paste0(pvals_dir, "/pvals_T200_laplace_delta", delta, "_h", h, "_s", s, ".csv"))[,-1]
pvals <- pvals[!is.na(rowSums(pvals)),]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:5 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

# Plot p-values against theoretical quantities from U(0,1)
dat3 <- data.frame(x=rep(seq(0, 1, length.out=nrow(pvals_sorted)), 4), y=as.vector(pvals_sorted[,c(1,2,3,5)]), 
                   N=rep(as.factor(c(1, 5, 10, 50)), each=nrow(pvals_sorted)))

fig14c <- ggplot(dat3) + geom_abline(intercept=0, slope=1) + 
  geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  labs(x="U(0,1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

###################################################################################################################
# Create Figure 15

pdf(paste0(plots_dir, "/fig14_laplace_qq.pdf"), width=20, height=6)
ggarrange(fig14a, fig14b, fig14c, ncol=3, common.legend=TRUE, legend="bottom")
dev.off()