#!/usr/bin/env Rscript

# generate Figure 4 based on results from optimize.py
source(file="03_figures/params.r")
fitness <- function(x, shape, scale, location, alpha)
{
	return (1-alpha*pgamma(x-location, shape=shape, scale=scale))
}

cuts <- seq(0, 50, 0.5)

pdf("03_figures/fitnessCurve.pdf", family="Times")
par(cex=1.5, mar=c(5,5,2,2)+0.1)
plot(fitness(cuts, shape, scale, location, alpha) ~ cuts, xlim=c(0,20), ylim=c(0,1), type='l', 
     xlab='False positive rate', ylab='Relative fitness')
abline(v=3.5, col='red', lty=2)
dev.off()
