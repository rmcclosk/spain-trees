#!/usr/bin/env Rscript

# Detect cross-contamination.

d <- read.table("01_annotate.dat", header=T)

# Counts by sequence, within each sample (patient/day).
sample.counts <- aggregate(count~patient+day+seq, d, sum)

# For each patient, maximum number of times the sequence appears
# in any of their samples.
max.counts <- aggregate(count~patient+seq, sample.counts, max)

# http://stackoverflow.com/questions/2453326/fastest-way-to-find-second-third-highest-lowest-value-in-vector-or-column
second.largest <- function (x) {
    n <- length(x)
    if (n == 1) 0
    else sort(x, partial=n-1)[n-1]
}

# The sequence's maximum count will be from its patient of origin. If it occurs
# in a second sample, that's cross-contamination.
max.count <- aggregate(count~seq, max.counts, max)
max.count.elsewhere <- aggregate(count~seq, max.counts, second.largest)

count.pairs <- merge(max.count, max.count.elsewhere, by=c("seq"), suffixes=c(".primary", ".secondary"))
count.pairs <- count.pairs[count.pairs$count.secondary > 0,]
write.csv(count.pairs[,c("count.primary", "count.secondary")], file="02_crosscontam.csv", row.names=F, quote=F)


# Plot a kernel density of total counts, and one of cross-contaminated counts.
xmax <- 20
cutoff <- 16
all.counts <- d[d$count <= xmax, "count"]
contam.counts <- max.count.elsewhere[max.count.elsewhere$count > 0, "count"]

pdf("02_crosscontam.pdf", family="Times")
hist(contam.counts, col="red", lwd=2, freq=F, xlim=c(0, xmax), ylim=c(0, 0.3),
     main="cutoff for cross-contamination filter", xlab="count", ylab="frequency")
lines(density(all.counts), col="blue", lwd=2)
abline(v=cutoff, lty=2)
legend("topright", legend=c("all", "cross-contaminated"), fill=c("blue", "red"),
       bg="white")
dev.off()

# Write a new table containing cross-contamination information.
d <- merge(d, max.count.elsewhere, by=c("seq"), suffixes=c("", ".contam"))
write.table(d, "02_crosscontam.dat", quote=F, row.names=F, col.names=T)
