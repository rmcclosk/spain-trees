#!/usr/bin/env Rscript

library(seqinr)
source(file="constants.r")

d <- read.table("01_annotate.dat", header=T)

# Get only the sequences from the first and last time point.
day.extremes <- rbind(aggregate(day~patient, d, max), aggregate(day~patient, d, min))
d <- merge(d, day.extremes)

# Remove sequences whose total count (in an individual sample) is less than the
# cutoff (min.count, defined in constants.r).
sample.counts <- aggregate(count~seq+patient+day, d, sum)
d <- merge(d, sample.counts, by=c("seq", "patient", "day"), suffixes=c("", ".sample"))
d <- d[d$count.sample >= min.count,] 

# Join the names of identical sequences with a "|", and remove duplicates.
names <- aggregate(taxa~seq+patient, d, paste, collapse="|")
d <- merge(d, names, by=c("seq", "patient"), suffixes=c("", ".collapsed"))

# Write a data file, to make it easier later.
stopifnot(duplicated(d$seq) == duplicated(cbind(d$seq, d$patient)))
stopifnot(duplicated(d$seq) == duplicated(d$taxa.collapsed))
write.table(d, "04_fasta.dat", col.names=T, row.names=F, quote=F)

# Write the sequences to fasta files, excluding duplicates.
d <- d[!duplicated(d$seq),]
d$seq <- as.character(d$seq)
. <- by(cbind(d$seq, d$taxa.collapsed, d$patient), d$patient, function (x) {
    fn <- sprintf("04_fasta/%02d.fasta", as.integer(as.character(x$V3[1])))
    seqs <- strsplit(as.character(x$V1), split=NULL)
    write.fasta(sequences=seqs, names=x$V2, file.out=fn)
})

