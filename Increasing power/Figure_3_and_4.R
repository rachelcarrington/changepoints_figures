# Figures 3 and 4

set.seed(100)

n <- 1000
num_iter <- 1000

for ( delta in c(1, 2, 3) ){
  for ( h in c(10, 20, 50) ){
    pvals <- matrix(NA, nrow=num_iter, ncol=6)
    iter <- 1
    while( iter <= num_iter ){
      x <- rnorm(n) + c(rep(0, n/2), rep(delta, n/2))
      results <- find_changepoints(x, "bs", list(threshold=3, maxiter=1))
      if ( length(results$changepoints) >= 1 ){
        if ( results$changepoints[1] >= h & results$changepoints[1] <= n - h ){
          pval_results <- calculate_pvals_all(results, N=50, h=h, sigma2=1, include_original=TRUE, num_pvals=1, return_probs=TRUE)
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
    }
    pvals[pvals > 1] <- 1 # some p-values are slightly above 1 because of floating point errors
    colnames(pvals) <- paste0("N", c(1, 2, 5, 10, 20, 50))
    write.csv(pvals, paste0("pvals_T", n, "_bs_delta", delta, "_h", h, ".csv"))
  }
}


###################################################################################################################
# Figure 3: QQ plots for BS under H1
###################################################################################################################
# Figure 3(a)

pvals <- read.csv("pvals_T1000_bs_delta1_h10.csv")[,-1]

pvals_sorted <- matrix(NA, nrow=1000, ncol=ncol(pvals))
for ( i in 1:ncol(pvals) ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

dat <- data.frame(x=rep(seq(0, 1, length.out=nrow(pvals_sorted)), ncol(pvals_sorted)), y=c(pvals_sorted), 
                  N=as.factor(rep(c(1, 2, 5, 10, 20, 50), each=nrow(pvals_sorted))))

fig3a <- ggplot(dat) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + 
  labs(x="U(0,1)", y="P-value", colour="Number of samples", linetype="Number of samples") + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22), legend.title=element_text(size=22), legend.text=element_text(size=22))

###################################################################################################################
# Figure 3(b)

pvals <- read.csv("pvals_T1000_bs_delta1_h10.csv")[,-1]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:6 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

dat <- data.frame(x=c(rep(pvals_sorted[,1], 5), NA), y=c(pvals_sorted[,2:6], NA), 
                  N=as.factor(c(rep(c(2, 5, 10, 20, 50), each=nrow(pvals_sorted)), 1)))

fig3b <- ggplot(dat) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="P-value (N = 1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

###################################################################################################################
# Figure 3(c)

pvals <- read.csv("pvals_T1000_bs_delta1_h50.csv")[,-1]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:ncol(pvals) ){
  pvals_sorted[,i] <- sort(pvals[,i])
}
dat <- data.frame(x=c(rep(pvals_sorted[,1], 5), NA), y=c(pvals_sorted[,2:6], NA), 
                  N=as.factor(c(rep(c(2, 5, 10, 20, 50), each=nrow(pvals_sorted)), 1)))

fig3c <- ggplot(dat) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="P-value (N = 1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))


###################################################################################################################
# Figure 3(d)

pvals <- read.csv("pvals_T1000_bs_delta2_h10.csv")[,-1]

pvals_sorted <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals))
for ( i in 1:ncol(pvals) ){
  pvals_sorted[,i] <- sort(pvals[,i])
}
dat <- data.frame(z=c(rep(pvals_sorted[,1], 5), NA), p=c(pvals_sorted[,2:6], NA), 
                  N=as.factor(c(rep(c(2, 5, 10, 20, 50), each=nrow(pvals_sorted)), 1)))

fig3d <- ggplot(dat) + geom_line(aes(x=z, y=p, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="P-value (N = 1)", y="") + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

###################################################################################################################
# Create Figure 3

ggarrange(fig3a + theme(legend.key.size=unit(4, "line")), fig3b, fig3c, fig3d, ncol=4, nrow=1, common.legend=TRUE, legend="bottom")

###################################################################################################################
# Figure 4: power plots for BS, K = 1
###################################################################################################################

N_list <- c(1, 2, 5, 10, 20)
power_results <- numeric(0)
for ( delta in c(1, 2, 3) ){
  for ( h in c(10, 20, 50) ){
    pvals <- read.csv(paste0("pvals_T1000_bs_delta", delta, "_h", h, ".csv"))[,-1]
    for ( i in 1:5 ){
      power_results <- rbind(power_results, c(h, delta, N_list[i], mean(pvals[,i] < 0.05)))
    }
  }
}
colnames(power_results) <- c("h", "delta", "N", "power")
power_results <- data.frame(power_results)
power_results$delta <- as.factor(power_results$delta)


fig4a <- ggplot(power_results[power_results$h==10,]) + geom_point(aes(x=N, y=power, colour=delta, shape=delta), size=3) + 
    geom_line(aes(x=N, y=power, colour=delta, linetype=delta), linewidth=2) +
    labs(x="Number of samples", y="Power") + coord_cartesian(ylim=c(0.25, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=22), axis.title=element_text(size=24))

fig4b <- ggplot(power_results[power_results$h==20,]) + geom_point(aes(x=N, y=power, colour=delta), size=3) + 
    geom_line(aes(x=N, y=power, colour=delta, linetype=delta), linewidth=2) +
    labs(x="Number of samples", y="") + coord_cartesian(ylim=c(0.25, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=22), axis.title=element_text(size=24))

fig4c <- ggplot(power_results[power_results$h==50,]) + geom_point(aes(x=N, y=power, colour=delta), size=3) + 
    geom_line(aes(x=N, y=power, colour=delta, linetype=delta), linewidth=2) +
    labs(x="Number of samples", y="") + coord_cartesian(ylim=c(0.25, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=22), axis.title=element_text(size=24))

fig4a <- fig4a + theme(legend.text=element_text(size=22), legend.title=element_text(size=22), legend.key.size=unit(3, "line")) +
         labs(colour="Size of change", shape="Size of change", linetype="Size of change")

pdf(paste0(plots_dir, "/fig4_power_plots_K1.pdf"), width=20, height=6)
ggarrange(fig4a, fig4b, fig4c, ncol=3, nrow=1, common.legend=TRUE, legend="bottom")
dev.off()