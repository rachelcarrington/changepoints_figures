library(changepoint)

# Get dataset
data("HC1")
n <- 2000
x <- HC1[1:n]

# Find changes and estimate variance
results <- find_changepoints(x, "bs", list(threshold=450))
sigma2 <- mad(x)

# Calculate p-values
set.seed(100)
pval_results_bs <- calculate_pvals_all(results, h=10, N=20, sigma2=sigma2, num_pvals=38, include_original=TRUE, return_probs=TRUE)
p1 <- pval_results_bs$P_both[,1] / pval_results_bs$P_phi_in_S[,1]
p10 <- rowSums(pval_results_bs$P_both[,1:10]) / rowSums(pval_results_bs$P_phi_in_S[,1:10])
p20 <- pval_results_bs$p_value

# Adjustment for multiple changes
p1_adjusted <- cbind(p1, p.adjust(p1, "bonferroni"), p.adjust(p1, "holm"), p.adjust(p1, "BH"))
p10_adjusted <- cbind(p10, p.adjust(p10, "bonferroni"), p.adjust(p10, "holm"), p.adjust(p10, "BH"))
p20_adjusted <- cbind(p20, p.adjust(p20, "bonferroni"), p.adjust(p20, "holm"), p.adjust(p20, "BH"))

###################################################################################################################
# Plots

g <- ggplot(data.frame(z=1:n, x=x)) + geom_point(aes(x=z, y=x)) 
cols1 <- cols10 <- cols20 <- rep("grey", 38)
cols1[p1_adjusted[,3] < 0.05] <- "red"
cols10[p10_adjusted[,3] < 0.05] <- "red"
cols20[p20_adjusted[,3] < 0.05] <- "red"

pdf(paste0(plots_dir, "/fig9_HC1_changes_bs_original_h10.pdf"), width=20, height=5)
g + geom_vline(xintercept = results$changepoints, colour=cols1) + labs(x="Position", y="GC content") + 
  theme_classic() + 
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))
dev.off()

pdf(paste0(plots_dir, "/fig9_HC1_changes_bs_new_h10_N10.pdf"), width=20, height=5)
g + geom_vline(xintercept = results$changepoints, colour=cols10) + labs(x="Position", y="GC content") + 
  theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))
dev.off()

###################################################################################################################
# Find out how many significant p-values there are

num_sig_pvals <- matrix( c(sum(p1_adjusted[,1] < 0.05), sum(p1_adjusted[,2] < 0.05), sum(p1_adjusted[,3] < 0.05), sum(p1_adjusted[,4] < 0.05),
  sum(p10_adjusted[,1] < 0.05), sum(p10_adjusted[,2] < 0.05), sum(p10_adjusted[,3] < 0.05), sum(p10_adjusted[,4] < 0.05),
  sum(p20_adjusted[,1] < 0.05), sum(p20_adjusted[,2] < 0.05), sum(p20_adjusted[,3] < 0.05), sum(p20_adjusted[,4] < 0.05)), nrow=3, byrow=TRUE)
rownames(num_sig_pvals) <- paste0("bs", c(1, 10, 20))
colnames(num_sig_pvals) <- c("raw", "bonferroni", "holm", "bh")
print(num_sig_pvals)
