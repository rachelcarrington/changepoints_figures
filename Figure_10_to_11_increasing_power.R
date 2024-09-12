# Increasing power

###########################################################################################################################################################################

# CUSUM: H0

set.seed(100)
h <- 10
n <- 200
NN <- 1000

p1 <- p5 <- p10 <- p20 <- p50 <- rep(NA, NN)
iter <- 1
while ( iter <= NN ){
  x <- rnorm(n)
  results <- find_changepoints(x, method="bs", model="var", params=list(threshold=1, maxiter=1, loss="cusum"))
  if ( length(results$results$b) >= 1 ){
    if ( results$results$b[1] >= h & results$results$b[1] <= n - h ){
      z <- calculate_pvals_all(results, h=h, N=50, sigma2=1, return_probs=TRUE)
      p1[iter] <- z$P_both[1] / z$P_phi_in_S[1]
      p5[iter] <- sum(z$P_both[1:5]) / sum(z$P_phi_in_S[1:5])
      p10[iter] <- sum(z$P_both[1:10]) / sum(z$P_phi_in_S[1:10])
      p20[iter] <- sum(z$P_both[1:20]) / sum(z$P_phi_in_S[1:20])
      p50[iter] <- sum(z$P_both[1:50]) / sum(z$P_phi_in_S[1:50])
      iter <- iter + 1
    }
  }
  print(iter)
}

p <- c(sort(p1), sort(p5), sort(p10), sort(p20), sort(p50))
saveRDS(p, "pvals_H0_psi_T200_h10.rds")

# Plot

g1 <- ggplot(data.frame(p=p, z=rep(seq(0, 1, length.out=NN), 5), NW=as.factor(rep(c(1, 5, 10, 20, 50), each=NN)))) +
  geom_line(aes(x=z, y=p, colour=NW, linetype=NW), linewidth=1.2) + 
  geom_abline(intercept=0, slope=1) + 
  theme_classic() + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22), legend.text=element_text(size=18), legend.title=element_text(size=18),
        legend.key.size=unit(3, "line")) +
  labs(x="U(0, 1)", y="p-value") + 
  guides(colour=guide_legend(title="Number of samples"), linetype=guide_legend(title="Number of samples"))

###########################################################################################################################################################################

# CUSUM: H1

set.seed(100)
h <- 10
n <- 200
NN <- 1000

p1 <- p5 <- p10 <- p20 <- p50 <- rep(NA, NN)
iter <- 1
while ( iter <= NN ){
  x <- c(rnorm(n/2), rnorm(n/2, sd=2))
  results <- find_changepoints(x, method="bs", model="var", params=list(threshold=1, maxiter=1, loss="cusum"))
  if ( length(results$results$b) >= 1 ){
    if ( results$results$b[1] >= h & results$results$b[1] <= n - h ){
      z <- calculate_pvals_all(results, h=h, N=50, sigma2=1, return_probs=TRUE)
      p1[iter] <- z$P_both[1] / z$P_phi_in_S[1]
      p5[iter] <- sum(z$P_both[1:5]) / sum(z$P_phi_in_S[1:5])
      p10[iter] <- sum(z$P_both[1:10]) / sum(z$P_phi_in_S[1:10])
      p20[iter] <- sum(z$P_both[1:20]) / sum(z$P_phi_in_S[1:20])
      p50[iter] <- sum(z$P_both[1:50]) / sum(z$P_phi_in_S[1:50])
      iter <- iter + 1
    }
  }
  print(iter)
}

p <- c(sort(p1), sort(p5), sort(p10), sort(p20), sort(p50))
saveRDS(p, "pvals_H1_cusum_psi_T200_h10.rds")

# Plot
g2 <- ggplot(data.frame(p=p, z=rep(seq(0, 1, length.out=NN), 5), N=as.factor(rep(c(1, 5, 10, 20, 50), each=NN)))) +
  geom_line(aes(x=z, y=p, colour=N, linetype=N), linewidth=1.2) + 
  geom_abline(intercept=0, slope=1) + 
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22), legend.text=element_text(size=18), legend.title=element_text(size=18),
        legend.key.size=unit(3, "line")) +
  labs(x="U(0, 1)", y="") + 
  guides(colour=guide_legend(title="Number of samples"), linetype=guide_legend(title="Number of samples"))

###########################################################################################################################################################################

# Figure 10: Plot CUSUM H0 and H1
ggarrange(g1, g2, common.legend=TRUE, legend="bottom")

###########################################################################################################################################################################

# Figure 11

# LRS: simulating under H0

set.seed(10)

NN <- 200
N_list <- c(50, 100, 200)
p1 <- p2 <- p5 <- p10 <- p20 <- matrix(NA, nrow=NN, ncol=length(N_list))
colnames(p1) <- colnames(p2) <- colnames(p5) <- colnames(p10) <- colnames(p20) <- N_list
h <- 20
n <- 200
N2 <- 20
threshold <- 2
maxiter <- 1

iter <- 1
while ( iter <= NN ){
  x <- rnorm(n)
  results <- find_changepoints(x, "bs", list(threshold=threshold, maxiter=maxiter, loss="lrs"), model="var")
  if ( length(results$results$b) >= 1 & results$results$b[1] >= h & results$results$b[1] <= length(x) - h ){
    b <- results$results$b[1]
    for ( k in 1:length(N_list) ){
      N <- N_list[k]
      pvals <- calculate_pvals_var(results, h, N_sample=N, NW=20)
      p1[iter, k] <- pvals$p_W[1]
      p2[iter, k] <- sum(pvals$p_W[1:2] * pvals$P_phi_in_S[1:2], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:2], na.rm=TRUE)
      p5[iter, k] <- sum(pvals$p_W[1:5] * pvals$P_phi_in_S[1:5], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:5], na.rm=TRUE)
      p10[iter, k] <- sum(pvals$p_W[1:10] * pvals$P_phi_in_S[1:10], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:10], na.rm=TRUE)
      p20[iter, k] <- sum(pvals$p_W[1:20] * pvals$P_phi_in_S[1:20], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:20], na.rm=TRUE)
    }
    iter <- iter + 1
    print(iter)
  }
}

all_pvalues <- list(p1=p1, p2=p2, p5=p5, p10=p10, p20=20)
saveRDS(all_pvalues, "pvalues_gp_lrs_H0_increasing_power.rds")

# Simulating under H1

set.seed(10)

NN <- 200
N_list <- c(50, 100, 200)
p1 <- p2 <- p5 <- p10 <- p20 <- matrix(NA, nrow=NN, ncol=length(N_list))
colnames(p1) <- colnames(p2) <- colnames(p5) <- colnames(p10) <- colnames(p20) <- N_list
h <- 20
n <- 200
N2 <- 20
threshold <- 2
maxiter <- 1

iter <- 1
while ( iter <= NN ){
  x <- c(rnorm(n/2), rnorm(n/2, sd=2))
  results <- find_changepoints(x, "bs", list(threshold=threshold, maxiter=maxiter, loss="lrs"), model="var")
  if ( length(results$results$b) >= 1 & results$results$b[1] >= h & results$results$b[1] <= length(x) - h ){
    b <- results$results$b[1]
    for ( k in 1:length(N_list) ){
      N <- N_list[k]
      pvals <- calculate_pvals_var(results, h, N_sample=N, NW=20)
      p1[iter, k] <- pvals$p_W[1]
      p2[iter, k] <- sum(pvals$p_W[1:2] * pvals$P_phi_in_S[1:2], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:2], na.rm=TRUE)
      p5[iter, k] <- sum(pvals$p_W[1:5] * pvals$P_phi_in_S[1:5], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:5], na.rm=TRUE)
      p10[iter, k] <- sum(pvals$p_W[1:10] * pvals$P_phi_in_S[1:10], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:10], na.rm=TRUE)
      p20[iter, k] <- sum(pvals$p_W[1:20] * pvals$P_phi_in_S[1:20], na.rm=TRUE) / sum(pvals$P_phi_in_S[1:20], na.rm=TRUE)
    }
    iter <- iter + 1
    print(iter)
  }
}

all_pvalues <- list(p1=p1, p2=p2, p5=p5, p10=p10, p20=p20)
saveRDS(all_pvalues, "pvalues_gp_lrs_H1_increasing_power.rds")

###########################################################################################################################################################################

# Figure 11: plots

pvals <- readRDS("pvalues_gp_lrs_H1_increasing_power.rds")
p1 <- pvals[[1]]
p2 <- pvals[[2]]
p5 <- pvals[[3]]
p10 <- pvals[[4]]
p20 <- pvals[[5]]

NN <- nrow(p1)

dat50 <- data.frame(p=c(sort(p1[,1]), sort(p2[,1]), sort(p5[,1]), sort(p10[,1]), sort(p20[,1])), 
                    z=rep(seq(0, 1, length.out=NN), 5), 
                    NW=as.factor(rep(c(1, 2, 5, 10, 20), each=NN)))

dat100 <- data.frame(p=c(sort(p1[,2]), sort(p2[,2]), sort(p5[,2]), sort(p10[,2]), sort(p20[,2])), 
                     z=rep(seq(0, 1, length.out=NN), 5), 
                     NW=as.factor(rep(c(1, 2, 5, 10, 20), each=NN)))

dat200 <- data.frame(p=c(sort(p1[,3]), sort(p2[,3]), sort(p5[,3]), sort(p10[,3]), sort(p20[,3])), 
                     z=rep(seq(0, 1, length.out=NN), 5), 
                     NW=as.factor(rep(c(1, 2, 5, 10, 20), each=NN)))


g50 <- ggplot(dat50) + 
  geom_line(aes(x=z, y=p, colour=NW, linetype=NW), linewidth=1.2) + 
  geom_abline(intercept=0, slope=1) + 
  theme_classic() + labs(x="U(0, 1)", y="p-value") +
  theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18), 
        legend.key.size=unit(3, "line")) + 
  guides(colour=guide_legend(title="Number of samples"), linetype=guide_legend(title="Number of samples"))

g100 <- ggplot(dat100) + 
  geom_line(aes(x=z, y=p, colour=NW, linetype=NW), linewidth=1.2) + 
  geom_abline(intercept=0, slope=1) + 
  theme_classic() + labs(x="U(0, 1)", y="") +
  theme(axis.title=element_text(size=22), axis.text=element_text(size=20)) + 
  guides(colour=guide_legend(title="Number of samples"), linetype=guide_legend(title="Number of samples"))

g200 <- ggplot(dat200) + 
  geom_line(aes(x=z, y=p, colour=NW, linetype=NW), linewidth=1.2) + 
  geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="U(0, 1)", y="") +
  theme(axis.title=element_text(size=22), axis.text=element_text(size=20)) + 
  guides(colour=guide_legend(title="Number of samples"), linetype=guide_legend(title="Number of samples"))


ggarrange(g50, g100, g200, ncol=3, common.legend=TRUE, legend="bottom")