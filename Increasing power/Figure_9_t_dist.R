# QQ plots, t distribution

###################################################################################################################
# H0

set.seed(1)

h_list <- c(10, 20)
df_list <- c(5, 10, 20)
n <- 200 # size of dataset
num_iter <- 1000 # number of simulations

for ( h in h_list ){
  for ( df in df_list ){
    set.seed(100)
    pvals <- matrix(NA, nrow=num_iter, ncol=6)
    iter <- 1
    while( iter <= num_iter ){
      x <- rt(n, df)
      results <- find_changepoints(x, "bs", list(threshold=3.5, maxiter=1))
      if ( length(results$changepoints) >= 1 ){ # only keep simulations for which the changepoint algorithm detects at least one changepoint
        if ( results$changepoints[1] >= h & results$changepoints[1] <= (n - h) ){
          results_pvals <- calculate_pvals_all(results, N=50, h=h, sigma2=df/(df - 2), num_pvals=1, include_original=TRUE, return_probs=TRUE)
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
    pvals[pvals > 1] <- 1 # some p-values may be slightly above 1 due to floating point errors
    pvals <- pvals[rowSums(is.na(pvals)) == 0,]
    colnames(pvals) <- paste0("N", c(1, 2, 5, 10, 20, 50))
    write.csv(pvals, paste0("pvals_tdist_H0_T", n, "_h", h, "_df", df, "_truevar.csv"))
  }
}

###################################################################################################################
# H1

set.seed(1)

num_iter <- 1000
n <- 200

for ( delta in c(1, 2) ){
  for ( h in c(10, 20) ){
    for ( df in c(5, 10, 20) ){
      pvals <- matrix(NA, nrow=num_iter, ncol=6)
      iter <- 1
      while( iter <= num_iter ){
        x <- rt(n, df) + c(rep(delta, n/2), rep(0, n/2))
        results <- find_changepoints(x, "bs", list(threshold=3.5, maxiter=1))
        if ( length(results$changepoints) >=1 ){
          if ( results$changepoints[1] >= h & results$changepoints[1] <= (n - h) ){
            results_pvals <- calculate_pvals_all(results, N=50, h=h, sigma2=df/(df - 2), include_original=TRUE, num_pvals=1, return_probs=TRUE)
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
      colnames(pvals) <- paste0("N", c(1, 2, 5, 10, 20, 50))
      write.csv(pvals, paste0("pvals_T200_tdist_delta", delta, "_h", h, "_df", df, ".csv"))
    }
  }
}

###################################################################################################################
# Fig. 8(a): H0, h = 10, df = 5

h <- 10
df <- 5
pvals <- read.csv(paste0(pvals_dir, "/pvals_tdist_H0_T200_h", h, "_df", df, "_truevar.csv"))[,-1]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:5 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

# Plot p-values against theoretical quantities from U(0,1)
dat <- data.frame(x=rep(seq(0, 1, length.out=nrow(pvals_sorted)), 4), y=as.vector(pvals_sorted[,c(1,2,3,5)]), 
                  N=rep(as.factor(c(1, 5, 10, 50)), each=nrow(pvals_sorted)))

fig8a <- ggplot(dat) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="P-value") + 
  theme_classic() +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=24), 
        legend.title=element_text(size=18), legend.text=element_text(size=18), legend.key.size=unit(3, "line"))

###################################################################################################################
# Fig. 8(b): H0, h = 10, df = 10

h <- 10
df <- 10
pvals <- read.csv(paste0(pvals_dir, "/pvals_tdist_H0_T200_h", h, "_df", df, "_truevar.csv"))[,-1]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:5 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

# Plot p-values against theoretical quantities from U(0,1)
dat <- data.frame(x=rep(seq(0, 1, length.out=nrow(pvals_sorted)), 4), y=as.vector(pvals_sorted[,c(1,2,3,5)]), 
                  N=rep(as.factor(c(1, 5, 10, 50)), each=nrow(pvals_sorted)))

fig8b <- ggplot(dat) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=24))

###################################################################################################################
# Fig. 8(c)

h <- 10
df <- 5
delta <- 1

pvals <- read.csv(paste0(pvals_dir, "/pvals_T200_tdist_delta", delta, "_h", h, "_df", df, ".csv"))[,-1]
pvals <- pvals[!is.na(rowSums(pvals)),]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:5 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

# Plot p-values against theoretical quantities from U(0,1)
dat <- data.frame(x=rep(seq(0, 1, length.out=nrow(pvals_sorted)), 4), y=as.vector(pvals_sorted[,c(1,2,3,5)]), 
                  N=rep(as.factor(c(1, 5, 10, 50)), each=nrow(pvals_sorted)))

fig8c <- ggplot(dat) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=24))

###################################################################################################################
# Create figure

pdf(paste0(plots_dir, "/fig8_t_dist_qq.pdf"), width=20, height=6)
ggarrange(fig8a, fig8b, fig8c, ncol=3, common.legend=TRUE, legend="bottom")
dev.off()