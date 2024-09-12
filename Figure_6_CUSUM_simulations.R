# Simulating under H0, fitting one change

h_list <- c(10, 20, 50)
n <- 200
N <- 1000

for ( h in h_list ){
  set.seed(1)
  p <- rep(NA, N)
  iter <- 1
  while ( iter <= N ){
    x <- rnorm(n)
    results <- find_changepoints(x, model="var", method="bs", params=list(loss="cusum", threshold=1, maxiter=1))
    if ( length(results$changepoints) >= 1 ){
      pval_results <- calculate_pvals_all(results, h, sigma2=1)
      p[iter] <- pval_results$p_value
      print(iter)
      iter <- iter + 1
    }
  }
  saveRDS(p, paste0("pvals_var_cusum_H0_T", n, "_h", h, ".rds"))
}


set.seed(1)
n <- 200
N <- 1000
p <- rep(NA, N)
iter <- 1
while ( iter <= N ){
  x <- rnorm(n)
  results <- find_changepoints(x, model="var", method="bs", params=list(loss="cusum", threshold=1, maxiter=1))
  if ( length(results$b) >= 1 ){
    pval_results <- calculate_pvals_all(results, h, sigma2=1)
    p[iter] <- pval_results$p_value
    print(iter)
    iter <- iter + 1
  }
}
saveRDS(p, paste0("pvals_var_cusum_H0_T", n, "_hNA.rds"))

###########################################################################################################################################################################

# Simulating under H1 with one change, fitting a single changepoint

set.seed(100)

sizes_of_change <- c(0.5, sqrt(2), 2)
h_list <- c(10, 20, 50, 0)
n <- 200 # or 1000
N <- 1000
p <- array(NA, dim=c(N, length(h_list), length(sizes_of_change)))
for ( k in 1:length(sizes_of_change) ){
  delta <- sizes_of_change[k]
  iter <- 1
  while ( iter <= N ){
    changepoint_loc <- sample(1:n, 1)
    x <- c(rnorm(changepoint_loc), rnorm(n - changepoint_loc, sd=delta))
    results <- find_changepoints(x, model="var", method="bs", params=list(loss="cusum", threshold=5, maxiter=1))
    if ( length(results$b) >= 1 ){
      for ( j in 1:length(h_list) ){
        if ( h_list[j] == 0 ){
          h <- NULL
        } else {
          h <- h_list[j]
        }
        pval_results <- calculate_pvals_all(results, h, sigma2=1, num_pvals=1)
        p[iter, j, k] <- pval_results$p_value
      }
      print(iter)
      iter <- iter + 1
    }
  }
}

dimnames(p)[[2]] <- h_list
dimnames(p)[[3]] <- sizes_of_change^2

saveRDS(p, paste0("pvals_var_cusum_H1_T", n, "_K1.rds"))

###########################################################################################################################################################################

# Plots

N <- 1000
n <- 200
p10 <- readRDS(paste0("pvals_var_cusum_H0_T", n, "_h10.rds"))
p20 <- readRDS(paste0("pvals_var_cusum_H0_T", n, "_h20.rds"))
p50 <- readRDS(paste0("pvals_var_cusum_H0_T", n, "_h50.rds"))
pNA <- readRDS(paste0("pvals_var_cusum_H0_T", n, "_hNA.rds"))

dat_H0 <- data.frame(z=rep(seq(0, 1, length.out=N)), p=c(sort(p10), sort(p20), sort(p50), sort(pNA)),
                     h=as.factor(rep(c("10", "20", "50", "NA"), each=N)))

p <- readRDS(paste0("pvals_var_cusum_H1_T", n, "_K1.rds"))
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
dat$h[dat$h==0] <- "NA"
dat$h <- as.factor(dat$h)

g0 <- ggplot(dat_H0) + 
  geom_abline(intercept=0, slope=1) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) + 
  theme_classic() +
  theme(legend.position="None", axis.title=element_text(size=28), axis.text=element_text(size=22), legend.title=element_text(size=26), 
        legend.text=element_text(size=24), legend.key.size=unit(3, "line")) +
  labs(x="U(0,1)", y="p-value")

g1 <- ggplot(dat[dat$delta == 1,]) + 
  geom_abline(intercept=0, slope=1) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) +
  theme_classic() + 
  theme(legend.position="None", axis.title=element_text(size=28), axis.text=element_text(size=22), legend.title=element_text(size=26), 
        legend.text=element_text(size=24), legend.key.size=unit(3, "line")) +
  labs(x="U(0,1)", y="")

g2 <- ggplot(dat[dat$delta == 2,]) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) +
  geom_abline(intercept=0, slope=1) + 
  theme_classic() +
  theme(legend.position="None", axis.title=element_text(size=28), axis.text=element_text(size=22), legend.title=element_text(size=26), 
        legend.text=element_text(size=24), legend.key.size=unit(3, "line")) +
  labs(x="U(0,1)", y="")

g3 <- ggplot(dat[dat$delta == 3,]) + 
  geom_line(aes(x=z, y=p, colour=h, linetype=h), linewidth=1.5) +
  geom_abline(intercept=0, slope=1) + 
  theme_classic() +
  theme(legend.position="None", axis.title=element_text(size=28), axis.text=element_text(size=22), legend.title=element_text(size=26), 
        legend.text=element_text(size=24), legend.key.size=unit(3, "line")) +
  labs(x="U(0,1)", y="")

ggarrange(g0, g1, g2, g3, ncol=4, common.legend=TRUE, legend="bottom")