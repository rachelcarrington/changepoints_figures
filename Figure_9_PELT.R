# Estimating changepoints using PELT 

library(ggplot2)
library(ggpubr)

#############################################################################################################################

## H0

set.seed(2)

NN <- 1000
N_list <- c(50, 100, 200)
fracs <- c(0.2, 0.5, 0.8, 1)
p_hat <- matrix(NA, nrow=NN, ncol=length(N_list) * length(fracs))
colnames(p_hat) <- paste0(rep(N_list, each=length(fracs)), "_", rep(c("gp", "both"), length(N_list)))
h <- 20
n <- 200

iter <- 1
while ( iter <= NN ){
  
  x <- rnorm(n)
  results <- find_changepoints(x, "pelt", list(penalty="Manual", pen.value=10), model="var")
  if ( length(results$changepoints) >= 1 & results$changepoints[1] >= h & results$changepoints[1] <= length(x) - h ){
    for ( k in 1:length(N_list) ){
      N <- N_list[k]
      
      for ( kk in 1:length(fracs) ){
        frac <- fracs[kk]
        p_hat[iter, (k - 1) * length(fracs) + kk] <- calculate_pvals_var(results, h, N, frac, num_pvals=1, cp_bound=FALSE)$p
      }
    }
    iter <- iter + 1
    print(iter)
  }
}

saveRDS(p_hat, paste0("pvalues_gp_lrs_pelt_H0_T", n, "_h", h, "_thresh10.rds"))

#############################################################################################################################

## H1 (with three equidistant changes) 

set.seed(2)

NN <- 1000
N_list <- c(50, 100, 200)
fracs <- c(0.2, 0.5, 0.8, 1)
p_hat <- matrix(NA, nrow=NN, ncol=length(N_list) * length(fracs))
colnames(p_hat) <- paste0(rep(N_list, each=length(fracs)), "_", rep(fracs, length(N_list)))
h <- 20
n <- 200
threshold <- 10

iter <- 1
while ( iter <= NN ){
  x <- c(rnorm(n/4), rnorm(n/4, sd=2), rnorm(n/4), rnorm(n/4, sd=0.5))
  results <- find_changepoints(x, "pelt", list(penalty="Manual", pen.value=threshold), model="var")
  if ( length(results$changepoints) >= 1 & results$changepoints[1] >= h & results$changepoints[1] <= length(x) - h ){
    for ( k in 1:length(N_list) ){
      N <- N_list[k]
      for ( kk in 1:length(fracs) ){
        frac <- fracs[kk]
        p_hat[iter, (k - 1) * length(fracs) + kk] <- calculate_pvals_var(results, h, N, frac, num_pvals=1, cp_bound=FALSE)$p
      }
    }
    iter <- iter + 1
    print(iter)
  }
}

saveRDS(p_hat, paste0("pvalues_gp_lrs_pelt_H1_T", n, "_h", h, "_thresh", threshold, "_K3_fixed_changes.rds"))

#############################################################################################################################

## H1 (3 random changes)

set.seed(2)

NN <- 1000
N_list <- c(50, 100, 200)
fracs <- c(0.2, 0.5, 0.8, 1)
p_hat <- matrix(NA, nrow=NN, ncol=length(N_list) * length(fracs))
colnames(p_hat) <- paste0(rep(N_list, each=length(fracs)), "_", rep(fracs, length(N_list)))
h <- 20
n <- 200
threshold <- 10

iter <- 1
while ( iter <= NN ){
  change_locs <- sort(sample(1:n, 3))
  while ( sum( c(change_locs, n) - c(0, change_locs) < h ) > 0 ){
    change_locs <- sort(sample(1:n, 3))
  }
  x <- c(rnorm(change_locs[1]), rnorm(change_locs[2] - change_locs[1], sd=2), rnorm(change_locs[3] - change_locs[2]), rnorm(n - change_locs[3], sd=0.5))
  results <- find_changepoints(x, "pelt", list(penalty="Manual", pen.value=threshold), model="var")
  if ( length(results$changepoints) >= 1 & results$changepoints[1] >= h & results$changepoints[1] <= length(x) - h ){
    for ( k in 1:length(N_list) ){
      N <- N_list[k]
      for ( kk in 1:length(fracs) ){
        frac <- fracs[kk]
        p_hat[iter, (k - 1) * length(fracs) + kk] <- calculate_pvals_var(results, h, N, frac, num_pvals=1, cp_bound=FALSE)$p
      }
    }
    iter <- iter + 1
    print(iter)
  }
}

saveRDS(p_hat, paste0("pvalues_gp_lrs_pelt_H1_T", n, "_h", h, "_thresh", threshold, "_K3_random_changes.rds"))

#############################################################################################################################

## Plots

p1 <- readRDS("pvalues_gp_lrs_pelt_H0_T200_h20_thresh10.rds")
fracs <- c(0.2, 0.5, 0.8, 1)
N_list <- c(50, 100, 200)
NN <- nrow(p1)
p_all <- z_all <- method_all <- N_all <- numeric(0)
for ( i in 1:ncol(p1) ){
  p_all <- c(p_all, sort(p1[,i]))
  z_all <- c(z_all, seq(0, 1, length.out=length(sort(p1[,i]))))
  method_all <- c(method_all, rep(fracs[ifelse(i %% length(fracs) == 0, length(fracs), i %% length(fracs))], length(sort(p1[,i]))))
  N_all <- c(N_all, rep(N_list[(i - 1) %/% length(fracs) + 1], length(sort(p1[,i]))))
}
dat1 <- data.frame(p=p_all, z=z_all, Method=as.factor(method_all), N=N_all)

g1 <- ggplot(dat1[dat1$N==50,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="", y="p-value", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))

g2 <- ggplot(dat1[dat1$N==100,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="", y="", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))

g3 <- ggplot(dat1[dat1$N==200,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="", y="", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))




p2 <- readRDS("pvalues_gp_lrs_pelt_H1_T200_h20_thresh10_K3_fixed_changes.rds")
fracs <- c(0.2, 0.5, 0.8, 1)
N_list <- c(50, 100, 200)
NN <- nrow(p2)
p_all <- z_all <- method_all <- N_all <- numeric(0)
for ( i in 1:ncol(p2) ){
  p_all <- c(p_all, sort(p2[,i]))
  z_all <- c(z_all, seq(0, 1, length.out=length(sort(p2[,i]))))
  method_all <- c(method_all, rep(fracs[ifelse(i %% length(fracs) == 0, length(fracs), i %% length(fracs))], length(sort(p2[,i]))))
  N_all <- c(N_all, rep(N_list[(i - 1) %/% length(fracs) + 1], length(sort(p2[,i]))))
}
dat2 <- data.frame(p=p_all, z=z_all, Method=as.factor(method_all), N=N_all)

g4 <- ggplot(dat2[dat2$N==50,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="", y="p-value", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))

g5 <- ggplot(dat2[dat2$N==100,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="", y="", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))

g6 <- ggplot(dat2[dat2$N==200,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="", y="", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))



p3 <- readRDS("pvalues_gp_lrs_pelt_H1_T200_h20_thresh10_K3_random_changes.rds")
fracs <- c(0.2, 0.5, 0.8, 1)
N_list <- c(50, 100, 200)
NN <- nrow(p3)
p_all <- z_all <- method_all <- N_all <- numeric(0)
for ( i in 1:ncol(p3) ){
  p_all <- c(p_all, sort(p3[,i]))
  z_all <- c(z_all, seq(0, 1, length.out=length(sort(p3[,i]))))
  method_all <- c(method_all, rep(fracs[ifelse(i %% length(fracs) == 0, length(fracs), i %% length(fracs))], length(sort(p3[,i]))))
  N_all <- c(N_all, rep(N_list[(i - 1) %/% length(fracs) + 1], length(sort(p3[,i]))))
}
dat3 <- data.frame(p=p_all, z=z_all, Method=as.factor(method_all), N=N_all)

g7 <- ggplot(dat3[dat3$N==50,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="U(0, 1)", y="p-value", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))

g8 <- ggplot(dat3[dat3$N==100,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="U(0, 1)", y="", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))

g9 <- ggplot(dat3[dat3$N==200,]) + geom_line(aes(x=z, y=p, colour=Method, linetype=Method), linewidth=1.2) +
    geom_abline(intercept=0, slope=1) + theme_classic() + labs(x="U(0, 1)", y="", colour="Proportion of samples", linetype="Proportion of samples") +
    theme(axis.title=element_text(size=22), axis.text=element_text(size=20), legend.title=element_text(size=18), legend.text=element_text(size=18))

ggarrange(g1 + theme(legend.key.size=unit(3, "line")), g2, g3, g4, g5, g6, g7, g8, g9, ncol=3, nrow=3, common.legend=TRUE, legend="bottom")
