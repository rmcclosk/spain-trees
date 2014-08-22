#!/usr/bin/env Rscript

# Reshape the CSV file to be more useful.

d <- read.csv("00_raw.csv", header=T)

d$day <- factor(sapply(as.character(d$Time.point), function (x) {
    strsplit(x, " ")[[1]][2]
}))
d$taxa <- with(d, paste(ENUM, ID, Direction, day, sep="_"))

keep.col <- c("taxa", "Patient", "day", "Count", "String.no.dashes",
              "G2P_FPR.")
d <- d[,keep.col]
colnames(d) <- c("taxa", "patient", "day", "count", "seq", "fpr")

# Remove patient 7 day 2, and also any sequences which g2p couldn't score.
d <- with(d, d[ (patient != 7 | day != 2) & !is.na(fpr) ,])

write.table(d, "01_annotate.dat", quote=F, col.names=T, row.names=F)
