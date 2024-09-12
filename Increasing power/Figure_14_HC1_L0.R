# Figure 13: GC content data plots, L0 segmentation

library(changepoint)
library(ChangepointInference)

data("HC1")
n <- 2000
x <- HC1[1:n]

x <- x / sd(x)
results <- changepoint_estimates(x, "L0", tuning_parameter=4.2)
sigma2 <- mad(x)

# 4.2 gives 38 changepoints
# 3.5 gives 46
# 2.7 gives 58
# 2 gives 77

set.seed(100)
results_l0 <- l0_segmentation_psi(y, lambda=4.2, N=20, h=10, sigma2=sigma2, sig=5, include_original=TRUE, num_pvals=38)
  # runs in about 16 minutes

p1 <- results_l0$p_value_orig
p10 <- rowSums(cbind(results_l0$P_both_orig, results_l0$P_both[,1:9])) / rowSums(cbind(results_l0$P_phi_in_S_orig, results_l0$P_phi_in_S[,1:9]))
p20 <- results_l0$p_value

p1_adjusted <- cbind(p1, p.adjust(p1, "bonferroni"), p.adjust(p1, "holm"), p.adjust(p1, "BH"))
p10_adjusted <- cbind(p10, p.adjust(p10, "bonferroni"), p.adjust(p10, "holm"), p.adjust(p10, "BH"))
p20_adjusted <- cbind(p20, p.adjust(p20, "bonferroni"), p.adjust(p20, "holm"), p.adjust(p20, "BH"))

g <- ggplot(data.frame(z=1:n, x=x)) + geom_point(aes(x=z, y=x)) + labs(x="", y="")
cols1 <- cols10 <- cols20 <- rep("grey", 38)
cols1[p1_adjusted[,3] < 0.05] <- "red"
cols10[p10_adjusted[,3] < 0.05] <- "red"
cols20[p20_adjusted[,3] < 0.05] <- "red"

pdf(paste0(plots_dir, "/fig13_HC1_changes_l0_original_h10.pdf"), width=20, height=5)
g + geom_vline(xintercept = results_l0$b, colour=cols1) + labs(x="Position", y="GC content") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))
dev.off()

pdf(paste0(plots_dir, "/fig13_HC1_changes_l0_new_h10_N10.pdf"), width=20, height=5)
g + geom_vline(xintercept = results_l0$b, colour=cols10) + labs(x="Position", y="GC content") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))
dev.off()

pdf(paste0(plots_dir, "/fig13_HC1_changes_l0_new_h10_N20.pdf"), width=20, height=5)
g + geom_vline(xintercept = results_l0$b, colour=cols20) + labs(x="Position", y="GC content") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))
dev.off()

# Find out how many significant p-values there are
num_sig_pvals <- matrix( c(sum(p1_adjusted[,1] < 0.05), sum(p1_adjusted[,2] < 0.05), sum(p1_adjusted[,3] < 0.05), sum(p1_adjusted[,4] < 0.05),
  sum(p10_adjusted[,1] < 0.05), sum(p10_adjusted[,2] < 0.05), sum(p10_adjusted[,3] < 0.05), sum(p10_adjusted[,4] < 0.05),
  sum(p20_adjusted[,1] < 0.05), sum(p20_adjusted[,2] < 0.05), sum(p20_adjusted[,3] < 0.05), sum(p20_adjusted[,4] < 0.05)), nrow=3, byrow=TRUE)
rownames(num_sig_pvals) <- paste0("l0_", c(1, 10, 20))
colnames(num_sig_pvals) <- c("raw", "bonferroni", "holm", "bh")
print(num_sig_pvals)
