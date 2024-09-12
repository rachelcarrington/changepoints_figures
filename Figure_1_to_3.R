library(ggplot2)

# Figure 1

set.seed(3)
h <- 30
x <- c(rnorm(40, sd=3), rnorm(60))
results <- binary_segmentation(x, model="var")
b <- results$changepoints[1]
phi_obs <- sum(x[(b - h + 1):b]^2) / sum(x[(b - h + 1):(b + h)]^2)
x2 <- calculate_x_phi_var(x, b, 0.5, phi_obs, h, h)

ggplot(data.frame(z=1:length(x), x=x, x2=x2)) + 
  geom_hline(yintercept=0, colour="grey") +
  geom_point(aes(x=z, y=x), size=2.5) +
  geom_point(aes(x=z, y=x2), colour="red", shape=1, size=2.5) +
  geom_vline(xintercept=b + 0.5, colour="red", linetype="dashed") +
  geom_vline(xintercept=c(b - h, b + h) + 0.5, colour="blue", linetype="dotted") +
  theme_classic() + labs(x="t", y="X") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

# Figure 2

set.seed(100)
x <- c(rnorm(200), rnorm(100, sd=sqrt(1.7)))
results <- find_changepoints(x, "bs", list(threshold=7, loss="cusum"), "var")

g1 <- ggplot(data.frame(z=1:300, x=x)) + 
  geom_line(aes(x=z, y=x)) + 
  geom_vline(xintercept=200, colour="blue") +
  geom_vline(xintercept=results$changepoints, colour="red", linetype="dashed") +
  theme_classic() + labs(x="", y="X") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

g2 <- ggplot(data.frame(z=1:300, x=x^2)) + 
  geom_line(aes(x=z, y=x)) + 
  geom_vline(xintercept=200, colour="blue") +
  geom_vline(xintercept=results$changepoints, colour="red", linetype="dashed") +
  theme_classic() + labs(x="", y="Y") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

g3 <- ggplot(data.frame(z=1:299, x=cusum(x^2))) + 
  geom_line(aes(x=z, y=x)) + 
  geom_vline(xintercept=200, colour="blue") +
  geom_vline(xintercept=results$changepoints, colour="red", linetype="dashed") +
  theme_classic() + labs(x="t", y="G(t)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22))

ggarrange(g1, g2, g3, nrow=3)

# Figure 3

set.seed(2)

NN <- 1000
cusum_cps <- lrs_cps <- matrix(NA, nrow=NN, ncol=3)
for ( iter in 1:NN ){
  x <- c(rnorm(100), rnorm(100, sd=2), rnorm(100, sd=0.5), rnorm(100))
  cusum_cps[iter,] <- sort(binary_segmentation(x, threshold=0, maxiter=3, model="var", loss="cusum")$changepoints)
  lrs_cps[iter,] <- sort(binary_segmentation(x, threshold=0, maxiter=3, model="var", loss="lrs")$changepoints)
}

props_found <- matrix(NA, nrow=3, ncol=2)
colnames(props_found) <- c("cusum", "lrs")
props_found[1, 1] <- mean(rowSums(abs(cusum_cps - 100) <= 10) >= 1)
props_found[2, 1] <- mean(rowSums(abs(cusum_cps - 200) <= 10) >= 1)
props_found[3, 1] <- mean(rowSums(abs(cusum_cps - 300) <= 10) >= 1)
props_found[1, 2] <- mean(rowSums(abs(lrs_cps - 100) <= 10) >= 1)
props_found[2, 2] <- mean(rowSums(abs(lrs_cps - 200) <= 10) >= 1)
props_found[3, 2] <- mean(rowSums(abs(lrs_cps - 300) <= 10) >= 1)
props_found

set.seed(2)
x <- c(rnorm(100), rnorm(100, sd=2), rnorm(100, sd=0.5), rnorm(100), rnorm(100, sd=0.5), rnorm(100))

results_bs <- binary_segmentation(x, threshold=6, maxiter=5, model="var", loss="cusum")
g1 <- ggplot(data.frame(z=1:length(x), x=x)) + 
  geom_point(aes(x=z, y=x)) + 
  geom_vline(xintercept=seq(100, length(x), by=100)) +
  geom_vline(xintercept=results_bs$changepoints, colour="red", linetype="dashed") +
  labs(x="", y="X") + theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

results_wbs <- wild_binary_segmentation(x, threshold=9, model="var", loss="cusum")
g2 <- ggplot(data.frame(z=1:length(x), x=x)) + 
  geom_point(aes(x=z, y=x)) + 
  geom_vline(xintercept=seq(100, length(x), by=100)) +
  geom_vline(xintercept=results_wbs$changepoints, colour="red", linetype="dashed") +
  labs(x="", y="X") + theme_classic() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20))

ggarrange(g1, g2, ncol=1)

