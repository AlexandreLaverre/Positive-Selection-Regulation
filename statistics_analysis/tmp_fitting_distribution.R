library(fitdistrplus)
library(actuar)
sp="human"
TF="FOXA1"

load(paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/tested_deltas_", sp, "_", TF, ".Rdata"))

mat <- as.matrix(deltas[["all.deltas.focal"]])
mat_obs <- as.matrix(deltas[["obs.deltas.focal"]])

mat.obs.data <- c(mat_obs[!is.na(mat_obs)])
mat.all.data <- c(mat[!is.na(mat)])
mat.obs.data.pos <- mat.obs.data+30
mat.all.data.pos <- sample(mat.all.data+30, 500000)
SVM.obs <- pval.higher$SVM+abs(round(min(pval.higher$SVM)))+1
deltaSVM.obs <- pval.higher$deltaSVM+abs(round(min(pval.higher$deltaSVM)))+1

descdist(mat.obs.data.pos, boot = 100)
descdist(mat.all.data.pos, boot = 100)

obs.fit <- list()
all.fit <- list()
SVM.fit <- list()
distributions <- c("norm", "gamma", "lnorm", "pareto", "weibull", "llogis", "burr")
start.values <- list("llogis" = list(shape = 1, scale = 500), "burr" = list(shape1 = 0.3, shape2 = 1, rate = 1))

for (dist in distributions){
  print(dist)
  if (dist == "burr"){
    obs.fit[[dist]] <- fitdist(mat.obs.data.pos, dist, start=start.values[[dist]], lower=c(0,0))
    all.fit[[dist]] <- fitdist(mat.all.data.pos, dist, start=start.values[[dist]], lower=c(0,0))
    SVM.fit[[dist]] <- fitdist(deltaSVM.obs, dist, start=start.values[[dist]], lower=c(0,0))
  }else{
    obs.fit[[dist]] <- fitdist(mat.obs.data.pos, dist, start=start.values[[dist]])
    all.fit[[dist]] <- fitdist(mat.all.data.pos, dist, start=start.values[[dist]])
    SVM.fit[[dist]] <- fitdist(deltaSVM.obs, dist, start=start.values[[dist]])
  }
  
  summary(obs.fit[[dist]])
  summary(all.fit[[dist]])
  summary(SVM.fit[[dist]])
}

pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/Fit_distributions_all_deltas.pdf")
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
denscomp(all.fit, legendtext = distributions, xlab="all substitutions")
qqcomp(all.fit, legendtext = distributions)
cdfcomp(all.fit, legendtext = distributions)
cdfcomp(all.fit, legendtext = distributions, xlogscale=T, ylogscale = T)
ppcomp(all.fit, legendtext = distributions)
gofstat(all.fit, fitnames = distributions)
dev.off()

# Observed deltas
pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/Fit_distributions_obs_deltas.pdf")
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
denscomp(obs.fit, legendtext = distributions, xlab="obs substitutions")
qqcomp(obs.fit, legendtext = distributions)
cdfcomp(obs.fit, legendtext = distributions)
cdfcomp(obs.fit, legendtext = distributions, xlogscale=T, ylogscale = T)
ppcomp(obs.fit, legendtext = distributions)
gofstat(obs.fit, fitnames = distributions)
dev.off()

### SVM
pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/Fit_distributions_SVM.pdf")
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
denscomp(SVM.fit, legendtext = distributions, xlab="obs SVM")
qqcomp(SVM.fit, legendtext = distributions)
cdfcomp(SVM.fit, legendtext = distributions)
cdfcomp(SVM.fit, legendtext = distributions, xlogscale=T, ylogscale = T)
ppcomp(SVM.fit, legendtext = distributions)
gofstat(SVM.fit, fitnames = distributions)
dev.off()


obs.delta <- mat.obs.data
all.delta <- sample(mat.all.data, 500000)

ggplot() + theme_light() +
  stat_qq(aes(sample = obs.delta), colour = "green") +
  stat_qq(aes(sample = all.delta), colour = "red") +
  stat_qq() +
  scale_color_manual(name='Q-Q plot',
                     breaks=c('all delta', 'obs delta'),
                     values=c('all delta'='red', 'obs delta'='green'))

gg_qq_empirical <- function(a, b, quantiles = seq(0, 1, 0.00001))
{
  a_lab <- deparse(substitute(a))
  if(missing(b)) {
    b <- rnorm(length(a), mean(a), sd(a))
    b_lab <- "normal distribution"
  }
  else b_lab <- deparse(substitute(b))
  
  ggplot(mapping = aes(x = quantile(a, quantiles), 
                       y = quantile(b, quantiles))) + 
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    labs(x = paste(deparse(substitute(a)), "quantiles"), 
         y = paste(deparse(substitute(b)), "quantiles"),
         title = paste("Empirical qq plot of", a_lab, "against", b_lab))
}

pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/all_versus_obs_deltas.pdf")
col=c(rgb(1,0,0,0.5), rgb(0,1,0,0.5))
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hist(mat, breaks=100, xlim=c(-20, 20), main=paste(sp, TF), freq=F, col=col[1], xlab="deltaSVM per mutation", ylim=c(0,0.20))
hist(mat_obs, breaks=100, xlim=c(-20, 20), col=col[2], freq=F, add=T)

legend("topright", legend=c("All sub", "Obs sub"), fill=col, bty="n")
abline(v=0, col="red")
mtext(paste("mean All=", signif(mean(mat, na.rm=T),3)), line=-10, at=15, cex=0.7)
mtext(paste("mean Obs=", signif(mean(mat_obs, na.rm=T),3)), line=-11, at=15, cex=0.7)

boxplot(list(mat, mat_obs), outline=F, notch=T, col=col, names=c("All sub", "Obs sub"), ylab="deltaSVM per mutation")

qq <- gg_qq_empirical(all.delta, obs.delta)
qq + theme_light() + coord_equal()
dev.off()

