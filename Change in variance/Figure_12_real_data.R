set.seed(12)
dat <- read.csv("features_target.csv")
x <- dat[,3]
N <- 100
h <- 50

# PELT

params <- list(penalty="Manual", pen.value=22)
results <- find_changepoints(x, "pelt", params, model="var")
p_hat <- calculate_pvals_var(results, h, N, cp_bound=FALSE)
p1 <- p.adjust(p_hat$p, method="holm")
b1 <- results$changepoints
g1 <- ggplot(data.frame(z=1:length(x), x=x)) + geom_point(aes(x=z, y=x)) +
  geom_vline(xintercept=b1, colour=ifelse(p1 < 0.05, "red", "grey")) +
  theme_classic() + labs(x="Index", y="X") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=22)) + 
  scale_x_continuous(expand = c(0.001,0)) + coord_cartesian(xlim=c(0, length(x)))

# CUSUM

results <- find_changepoints(x, method="bs", params=list(loss="cusum", threshold=0.0005, maxiter=11), model="var")
p_hat <- calculate_pvals_all(results, h, sigma2=1)
p2 <- p.adjust(p_hat$p_value, method="holm")
b2 <- results$changepoints
g2 <- ggplot(data.frame(z=1:length(x), x=x)) + geom_point(aes(x=z, y=x)) +
  geom_vline(xintercept=b2, colour=ifelse(p2 < 0.05, "red", "grey")) +
  theme_classic() + labs(x="", y="X") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=22)) + 
  scale_x_continuous(expand = c(0.001,0)) + coord_cartesian(xlim=c(0, length(x)))

# LRS BS

params <- list(loss="lrs", threshold=18, maxiter=11)
results <- find_changepoints(x, method="bs", params=params, model="var")
p_hat <- calculate_pvals_var(results, h, N, cp_bound=FALSE)
p3 <- p.adjust(p_hat$p, method="holm")
b3 <- results$changepoints
g3 <- ggplot(data.frame(z=1:length(x), x=x)) + geom_point(aes(x=z, y=x)) +
  geom_vline(xintercept=b3, colour=ifelse(p3 < 0.05, "red", "grey")) +
  theme_classic() + labs(x="Index", y="X") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=22)) + 
  scale_x_continuous(expand = c(0.001,0)) + coord_cartesian(xlim=c(0, length(x)))

# Create figure
ggarrange(g2, g3, g1, ncol=1, nrow=3)
