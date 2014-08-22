#!/usr/bin/env Rscript

# Generate frequency through time plots.

d <- read.table("02_crosscontam.dat", header=T)

# Remove cross-contaminants.
d <- d[d$patient != 7 | d$day != 2,]
d <- d[d$count.contam == 0,]

# Bin data by FPR.
d$bin <- cut(d$fpr, c(0, 2, 3.5, 10, 20, 100), include.lowest=T)

# Aesthetics.
col <- c("red", "orange", "forestgreen", "blue", "black", "magenta")
ylab <- "Frequency"
xlab <- "Days"
mar <- c(3,4,1,0)+0.2
mgp <- c(1.7,0.5,0)
cex <- 0.8

# Make frequency by day plots.
pdf("03_figures/frequencyByDay.pdf", family="Times")
main <- sapply(0:2, function(j) {
	par(mfrow=c(3, 3), mar=mar, mgp=mgp, cex=cex)
	sapply((j*9+1):(j*9+9), function (p) {
		pat.d <- d[d$patient == p,]
		t <- prop.table(table(pat.d$bin, pat.d$day), 2)
		plot(NULL, xlim=c(0, 7), ylim=c(0, 1), ylab=ylab, xlab=xlab,
			main=paste0("P", p))
		sapply(1:nrow(t), function (i) {
			lines(colnames(t), t[i,], col=col[i], pch=i, type="b")
		})
	})
})
dev.off()

# Make legend.
pdf("03_figures/legend.pdf", family="Times")
plot.new()
legend(x=0, y=1, legend=c(
	expression(0 <= ~paste(FPR <= 2)),
	expression(2 < ~paste(FPR <= 3.5)),
	expression(3.5 < ~paste(FPR <= 10)), 
	expression(10 < ~paste(FPR <= 20)),
	expression(20 < FPR)), pch=1:5, col=col, lty=1)
warnings()
dev.off()

# Make detailed frequency by day plots for particular patients.
d <- d[d$patient %in% c(7, 13, 17, 19),]
d <- d[d$fpr > 3 & d$fpr <= 15,]
d$bin <- cut(d$fpr, c(3, 5, 7, 9, 11, 13, 15), include.lowest=T)

pdf("03_figures/frequencyByDay_detailed.pdf", family="Times")
par(mfrow=c(2, 2), mar=mar, mgp=mgp, cex=cex)
. <- sapply(unique(d$patient), function (p) {
	pat.d <- d[d$patient == p,]
	t <- prop.table(table(pat.d$bin, pat.d$day), 2)
	plot(NULL, xlim=c(0, 7), ylim=c(0, 1), ylab=ylab, xlab=xlab,
		main=paste0("P", p))
	sapply(1:nrow(t), function (i) {
		lines(colnames(t), t[i,], col=col[i], pch=i, type="b")
	})
})
dev.off()

# Make legend.
pdf("03_figures/legend_detailed.pdf", family="Times")
plot.new()
legend(x=0, y=1, legend=c(
	expression(3 <= ~paste(FPR <= 5)),
	expression(5 < ~paste(FPR <= 7)),
	expression(7 < ~paste(FPR <= 9)), 
	expression(9 < ~paste(FPR <= 11)),
	expression(11 < ~paste(FPR <= 13)),
	expression(13 < ~paste(FPR <= 15))), pch=1:6, col=col, lty=1)
warnings()
dev.off()
