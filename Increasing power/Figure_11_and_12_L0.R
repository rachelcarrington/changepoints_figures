# Figure 11: QQ plots for L0 segmentation

library(ChangepointInference)

#########################################################################################
# Simulating under H1 with one changepoint

n <- 1000
delta_list <- c(1, 2, 3)
h_list <- c(10, 20, 30)
num_iter <- 1000

for ( delta in delta_list ){
  for ( h in h_list ){
    pvals <- matrix(NA, nrow=num_iter, ncol=6)
    iter <- 1
    while( iter <= num_iter ){
      x <- rnorm(n) + c(rep(delta/2, n/2), rep(-delta/2, n/2))
      fit <- changepoint_estimates(x, "L0", 10)
      if ( length(fit$change_pts) >= 1 ){
        if ( fit$change_pts[1] >= h & fit$change_pts[1] <= (n - h) ){
          pval_results <- l0_segmentation_psi(x, lambda=10, N=50, h=h, sigma2=1, include_original=TRUE, num_pvals=1, sig=20)
          pvals[iter,1] <- pval_results$p_value_orig
          pvals[iter,2] <- (sum(pval_results$P_both[1]) + pval_results$P_both_orig) / (sum(pval_results$P_phi_in_S[1]) + pval_results$P_phi_in_S_orig)
          pvals[iter,3] <- (sum(pval_results$P_both[1:4]) + pval_results$P_both_orig) / (sum(pval_results$P_phi_in_S[1:4]) + pval_results$P_phi_in_S_orig)
          pvals[iter,4] <- (sum(pval_results$P_both[1:9]) + pval_results$P_both_orig) / (sum(pval_results$P_phi_in_S[1:9]) + pval_results$P_phi_in_S_orig)
          pvals[iter,5] <- (sum(pval_results$P_both[1:19]) + pval_results$P_both_orig) / (sum(pval_results$P_phi_in_S[1:19]) + pval_results$P_phi_in_S_orig)
          pvals[iter,6] <- (sum(pval_results$P_both) + pval_results$P_both_orig) / (sum(pval_results$P_phi_in_S) + pval_results$P_phi_in_S_orig)
          iter <- iter + 1
          print(iter)
        }
      }
    }
    saveRDS(pvals, paste0(pvals_dir, "/pvals_T", n, "_l0_delta", delta, "_h", h, ".rds"))
  }
}

###################################################################################################################
# Plots

# Fig. 11(a): p-values against U(0, 1) for h = 10, delta = 1

pvals <- readRDS(paste0(pvals_dir, "/pvals_T1000_l0_delta1_h10.rds"))
pvals_sorted <- matrix(NA, nrow=1000, ncol=6)
for ( i in 1:6 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}
dat <- data.frame(x=rep(seq(0, 1, length.out=1000), 6), y=as.vector(pvals_sorted), 
                  N=rep(as.factor(c(1, 2, 5, 10, 20, 50)), each=1000))

fig10a <- ggplot(dat) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="U(0,1)", y="P-value", colour="Number of samples", linetype="Number of samples") +
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22),
  legend.title=element_text(size=22), legend.text=element_text(size=22), legend.key.size=unit(4, "line"), legend.position="bottom")

# Figure 11(b)

dat2 <- data.frame(x=c(rep(pvals_sorted[,1], 5), NA), y=c(as.vector(pvals_sorted[,2:6]), NA), 
                  N=as.factor(c(rep(c(2, 5, 10, 20, 50), each=1000), 1)))

fig10b <- ggplot(dat2) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="P-value (N = 1)", y="") +
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22))

# Figure 11(c)

pvals <- readRDS(paste0(pvals_dir, "/pvals_T1000_l0_delta1_h30.rds"))

pvals_sorted <- matrix(NA, nrow=1000, ncol=6)
for ( i in 1:6 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

dat2 <- data.frame(x=c(rep(pvals_sorted[,1], 5), NA), y=c(as.vector(pvals_sorted[,2:6]), NA), 
                  N=as.factor(c(rep(c(2, 5, 10, 20, 50), each=1000), 1)))

fig10c <- ggplot(dat2) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="P-value (N = 1)", y="") +
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22))

# Figure 11(d)

pvals <- readRDS(paste0(pvals_dir, "/pvals_T1000_l0_delta2_h10.rds"))

pvals_sorted <- matrix(NA, nrow=1000, ncol=6)
for ( i in 1:6 ){
  pvals_sorted[,i] <- sort(pvals[,i])
}

dat2 <- data.frame(x=c(rep(pvals_sorted[,1], 5), NA), y=c(as.vector(pvals_sorted[,2:6]), NA), 
                  N=as.factor(c(rep(c(2, 5, 10, 20, 50), each=1000), 1)))

fig10d <- ggplot(dat2) + geom_line(aes(x=x, y=y, colour=N, linetype=N), linewidth=1.2) +
  geom_abline(intercept=0, slope=1) + labs(x="P-value (N = 1)", y="") +
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22))

###################################################################################################################
# Create figure

pdf(paste0(plots_dir, "/fig10_qq_plots_L0.pdf"), width=20, height=6)
ggarrange(fig10a, fig10b, fig10c, fig10d, ncol=4, common.legend=TRUE, legend="bottom")
dev.off()

###################################################################################################################
# Figure 12: Power, L0 segmentation

N_list <- c(1, 2, 5, 10, 20)
power_results <- numeric(0)
for ( delta in c(1, 2, 3) ){
  for ( h in c(10, 20, 30) ){
    pvals <- readRDS(paste0(pvals_dir, "/pvals_T1000_l0_delta", delta, "_h", h, ".rds"))
    for ( i in 1:5 ){
      power_results <- rbind(power_results,
                             c(h, delta, N_list[i], sum(pvals[,i] < 0.1), sum(pvals[,i] < 0.05), sum(pvals[,i] < 0.01)))
    }
  }
}
colnames(power_results) <- c("h", "delta", "N", ".1", ".05", ".01")
power_results <- data.frame(power_results)
power_results$delta <- as.factor(power_results$delta)

fig11a <- ggplot(power_results[power_results$h==10,]) + geom_point(aes(x=N, y=X.05/1000, colour=delta, shape=delta), size=3) + 
    geom_line(aes(x=N, y=X.05/1000, colour=delta, linetype=delta), linewidth=2) +
    labs(x="Number of samples", y="Power", colour="Size of change", linetype="Size of change", shape="Size of change") + 
    coord_cartesian(ylim=c(0.25, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=22), 
          legend.text=element_text(size=20), legend.title=element_text(size=20), legend.key.size=unit(4, "line"))

fig11b <- ggplot(power_results[power_results$h==20,]) + geom_point(aes(x=N, y=X.05/1000, colour=delta, shape=delta), size=3) + 
    geom_line(aes(x=N, y=X.05/1000, colour=delta, linetype=delta), linewidth=2) +
    labs(x="Number of samples", y="") + coord_cartesian(ylim=c(0.25, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

fig11c <- ggplot(power_results[power_results$h==30,]) + geom_point(aes(x=N, y=X.05/1000, colour=delta, shape=delta), size=3) + 
    geom_line(aes(x=N, y=X.05/1000, colour=delta, linetype=delta), linewidth=2) +
    labs(x="Number of samples", y="") + coord_cartesian(ylim=c(0.25, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

pdf(paste0(plots_dir, "/fig11_l0_power.pdf"), width=20, height=6)
ggarrange(fig11a, fig11b, fig11c, ncol=3, common.legend=TRUE, legend="bottom")
dev.off()