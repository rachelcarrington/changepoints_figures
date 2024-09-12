# Simulating data with 3 changes, estimating changepoints using CUSUM

library(ggplot2)
library(ggpubr)

# Fixed changes

set.seed(100)

sizes_of_change <- c(sqrt(2), 2)
h_list <- c(10, 20)
n <- 200
K <- 4
threshold <- 5
N <- 1000
pvals <- numeric(0)
for ( k in 1:length(sizes_of_change) ){
  delta <- sizes_of_change[k]
  iter <- 1
  while ( iter <= N ){
    x <- rnorm(n) * c(rep(1, n/5), rep(delta, n/5), rep(1, n/5), rep(delta, n/5), rep(1, n/5))
    results <- find_changepoints(x, model="var", method="bs", params=list(loss="cusum", threshold=threshold, maxiter=K))
    if ( length(results$b) >= 1 ){
      for ( j in 1:length(h_list) ){
        h <- h_list[j]
        pval_results <- calculate_pvals_all(results, h=h, sigma2=1)
        if ( length(pval_results$p_value) == K ){
          pvals <- rbind(pvals, c(delta, h, pval_results$p_value))
        } else {
          pvals <- rbind(pvals, c(delta, h, pval_results$p_value, rep(NA, K - length(pval_results$p_value))))
        }
      }
      print(iter)
      iter <- iter + 1
    }
  }
}
saveRDS(pvals, "pvals_cusum_multiple_pvals_fixed_changes.rds")


# Random changes

set.seed(100)

sizes_of_change <- c(sqrt(2), 2)
h_list <- c(10, 20)
n <- 200
K <- 4
threshold <- 3
N <- 1000
pvals <- numeric(0)
for ( k in 1:length(sizes_of_change) ){
  delta <- sizes_of_change[k]
  for ( j in 1:length(h_list) ){
    h <- h_list[j]
    iter <- 1
    while ( iter <= N ){
      change_locs <- sort(sample(1:n, K))
      while ( sum( c(change_locs, n) - c(0, change_locs) < h ) > 0 ){
        change_locs <- sort(sample(1:n, K))
      }
      x <- rnorm(n) * c(rep(1, change_locs[1]), rep(delta, change_locs[2] - change_locs[1]), rep(1, change_locs[3] - change_locs[2]), rep(delta, change_locs[4] - change_locs[3]), rep(1, n - change_locs[4]))
      results <- find_changepoints(x, model="var", method="bs", params=list(loss="cusum", threshold=threshold, maxiter=K))
      if ( length(results$b) >= 1 ){
        pval_results <- calculate_pvals_all(results, h=h, sigma2=1)
        if ( length(pval_results$p_value) == K ){
          pvals <- rbind(pvals, c(delta, h, pval_results$p_value))
        } else {
          pvals <- rbind(pvals, c(delta, h, pval_results$p_value, rep(NA, K - length(pval_results$p_value))))
        }
      }
      print(iter)
      iter <- iter + 1
    }
  }
}

saveRDS(pvals, "pvals_cusum_multiple_pvals_fixed_changes.rds")

###########################################################################################################################################################################

# Plots

p <- readRDS("pvals_var_cusum_H1_T200_K4_random_changes_threshold5.rds")
p[p > 1] <- 1
h_list <- dimnames(p)[[2]]
sizes_of_change <- sqrt(as.numeric(dimnames(p)[[3]]))
p_all <- z_all <- h_all <- delta_all <- numeric(0)
for ( j in 1:length(h_list) ){
  for ( k in 1:length(sizes_of_change) ){
    p_all <- c(p_all, sort(p[,j,k]))
    z_all <- c(z_all, seq(0, 1, length.out=sum(!is.na(p[,j,k]))))
    h_all <- c(h_all, rep(h_list[j], sum(!is.na(p[,j,k]))))
    delta_all <- c(delta_all, rep(k, sum(!is.na(p[,j,k]))))
  }
}
dat <- data.frame(z=z_all, p=p_all, h=h_all, delta=delta_all)
dat$h <- as.factor(dat$h)

g1r <- ggplot(dat[dat$delta == 1,]) + 
  geom_abline(intercept=0, slope=1) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) +
  theme_classic() + 
  theme(axis.title=element_text(size=28), axis.text=element_text(size=22)) +
  labs(x="U(0,1)", y="")

g2r <- ggplot(dat[dat$delta == 2,]) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) +
  geom_abline(intercept=0, slope=1) + 
  theme_classic() + 
  theme(axis.title=element_text(size=28), axis.text=element_text(size=22)) +
  labs(x="U(0,1)", y="")

p <- readRDS("pvals_var_cusum_H1_T200_K4_fixed_changes_threshold5.rds")
p[p > 1] <- 1
h_list <- dimnames(p)[[2]]
sizes_of_change <- sqrt(as.numeric(dimnames(p)[[3]]))
p_all <- z_all <- h_all <- delta_all <- numeric(0)
for ( j in 1:length(h_list) ){
  for ( k in 1:length(sizes_of_change) ){
    p_all <- c(p_all, sort(p[,j,k]))
    z_all <- c(z_all, seq(0, 1, length.out=sum(!is.na(p[,j,k]))))
    h_all <- c(h_all, rep(h_list[j], sum(!is.na(p[,j,k]))))
    delta_all <- c(delta_all, rep(k, sum(!is.na(p[,j,k]))))
  }
}
dat <- data.frame(z=z_all, p=p_all, h=h_all, delta=delta_all)
dat <- dat[dat$h != "0" & dat$h != "50",]
dat$h <- as.factor(dat$h)

g1f <- ggplot(dat[dat$delta == 1,]) + 
  geom_abline(intercept=0, slope=1) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) +
  theme_classic() + 
  theme(axis.title=element_text(size=28), axis.text=element_text(size=22)) +
  labs(x="U(0,1)", y="p-value")

g2f <- ggplot(dat[dat$delta == 2,]) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) +
  geom_abline(intercept=0, slope=1) + 
  theme_classic() + 
  theme(axis.title=element_text(size=28), axis.text=element_text(size=22)) +
  labs(x="U(0,1)", y="")

ggarrange(g1f + theme(legend.title=element_text(size=26), legend.text=element_text(size=24), legend.key.size=unit(3, "line")), 
          g2f, g1r, g2r, ncol=4, common.legend=TRUE, legend="bottom")
