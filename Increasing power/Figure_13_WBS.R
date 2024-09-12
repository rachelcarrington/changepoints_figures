# Figure 13: power plots for WBS, K = 4

###################################################################################################################
# Simulations

set.seed(100)

num_iter <- 100

for ( delta in c(1, 2, 3) ){
  for ( h in c(10, 20, 50) ){
    pvals <- array(NA, dim=c(num_iter, 6, 4))
    iter <- 1
    while( iter <= num_iter ){
      x <- rnorm(1000) + c(rep(delta/2, 100), rep(-delta/2, 300), rep(delta/2, 100), rep(-delta/2, 200), rep(delta/2, 300))
      # changes at t = 100, 400, 500, 700
      results <- find_changepoints(x, "wbs", list(threshold=4, num_rand_samples=100))
      if ( length(results$changepoints) >= 1 ){
        if ( results$changepoints[1] >= h & results$changepoints[1] <= (1000 - h) ){
          pval_results <- calculate_pvals_all(results, N=30, h=h, include_original=TRUE, num_pvals=4, sigma2=1, return_probs=TRUE)
          ncp <- length(pval_results$p_value)
          pvals[iter, 1, ] <- c(pval_results$P_both[,1] / pval_results$P_phi_in_S[,1], rep(NA, 4 - ncp))
          pvals[iter, 2, ] <- c(rowSums(pval_results$P_both[,1:2]) / rowSums(pval_results$P_phi_in_S[,1:2]), rep(NA, 4 - ncp))
          pvals[iter, 3, ] <- c(rowSums(pval_results$P_both[,1:5]) / rowSums(pval_results$P_phi_in_S[,1:5]), rep(NA, 4 - ncp))
          pvals[iter, 4, ] <- c(rowSums(pval_results$P_both[,1:10]) / rowSums(pval_results$P_phi_in_S[,1:10]), rep(NA, 4 - ncp))
          pvals[iter, 5, ] <- c(rowSums(pval_results$P_both[,1:20]) / rowSums(pval_results$P_phi_in_S[,1:20]), rep(NA, 4 - ncp))
          pvals[iter, 6, ] <- c(pval_results$p_value, rep(NA, 4 - ncp))
          iter <- iter + 1
          print(iter)
        }
      }
    }
    pvals[pvals > 1] <- 1 # some p-values are slightly above 1 because of floating point errors
    colnames(pvals) <- paste0("N", c(1, 2, 5, 10, 20, 30))
    
    saveRDS(pvals, paste0(pvals_dir, "/pvals_T1000_wbs_K4_delta", delta, "_h", h, ".rds"))
  }
}

###################################################################################################################
# Plots

N_list <- c(1, 2, 5, 10, 20, 50)
power_results <- numeric(0)
for ( delta in c(1, 2, 3) ){
  for ( h in c(10, 20, 50) ){
    pvals <- readRDS(paste0(pvals_dir, "/pvals_T1000_wbs_K4_delta", delta, "_h", h, ".rds"))
    for ( i in 1:6 ){
      for ( j in 1:4 ){
        power_results <- rbind(power_results,
                               c( h, delta, N_list[i], j, 
                                  sum(pvals[,i,j] < 0.1, na.rm=TRUE) / sum(!is.na(pvals[,i,j])), 
                                  sum(pvals[,i,j] < 0.05, na.rm=TRUE) / sum(!is.na(pvals[,i,j])), 
                                  sum(pvals[,i,j] < 0.01, na.rm=TRUE) / sum(!is.na(pvals[,i,j])) ))
      }
    }
  }
}
colnames(power_results) <- c("h", "delta", "N", "Changepoint.index", ".1", ".05", ".01")
power_results <- data.frame(power_results)
power_results$delta <- as.factor(power_results$delta)
power_results$Changepoint.index <- as.factor(power_results$Changepoint.index)
power_results$grp <- as.factor(paste0(power_results$delta, power_results$Changepoint.index))

fig12a <- ggplot(power_results[power_results$h==10,]) + geom_point(aes(x=N, y=X.05, colour=delta, shape=Changepoint.index), size=3) + 
    geom_line(aes(x=N, y=X.05, colour=delta, group=grp), linewidth=1.2) +
    labs(x="Number of samples", y="Power", colour="Size of change", shape="Changepoint index") + 
    coord_cartesian(ylim=c(0.2, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=22), legend.text=element_text(size=20), legend.title=element_text(size=20),
          legend.key.size=unit(3, "line"))

fig12b <- ggplot(power_results[power_results$h==20,]) + geom_point(aes(x=N, y=X.05, colour=delta, shape=Changepoint.index), size=3) + 
    geom_line(aes(x=N, y=X.05, colour=delta, group=grp), linewidth=1.2) +
    labs(x="Number of samples", y="") + coord_cartesian(ylim=c(0.2, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

fig12c <- ggplot(power_results[power_results$h==50,]) + geom_point(aes(x=N, y=X.05, colour=delta, shape=Changepoint.index), size=3) + 
    geom_line(aes(x=N, y=X.05, colour=delta, group=grp), linewidth=1.2) +
    labs(x="Number of samples", y="") + coord_cartesian(ylim=c(0.2, 1), xlim=c(0, 20)) +
    theme_classic() + scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

pdf(paste0(plots_dir, "/fig12_power_plots_WBS_K4.pdf"), width=20, height=6)
ggarrange(fig12a, fig12b, fig12c, ncol=3, nrow=1, common.legend=TRUE, legend="bottom")
dev.off()