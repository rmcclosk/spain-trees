# Draw a tree for each patient.

library(ape)
source(file="constants.r")

d <- read.table("04_fasta.dat", header=T)
d$tropism <- ifelse(d$fpr <= x4.fpr, "X4", "R5")

# We will root the trees at the most prevalent baseline sequence.
roots <- with(subset(d, day==0), by(cbind(taxa.collapsed, count.sample), patient, function(x) {
    x[which.max(x$count.sample),"taxa.collapsed"]
}))
roots <- levels(d$taxa.collapsed)[c(roots)]

# Find which samples occured at which time points.
d <- merge(d, aggregate(day~patient, d, min), by=c("patient"), suffixes=c("", ".first.patient"))
d <- merge(d, aggregate(day~patient, d, max), by=c("patient"), suffixes=c("", ".last.patient"))
d <- merge(d, aggregate(day~taxa.collapsed, d, min), by=c("taxa.collapsed"), suffixes=c("", ".min.taxa"))
d <- merge(d, aggregate(day~taxa.collapsed, d, max), by=c("taxa.collapsed"), suffixes=c("", ".max.taxa"))

# Labels for first and last time points.
first <- subset(d, day==day.first.patient, select=c(taxa.collapsed, count.sample))
colnames(first)[2] <- "label.first"
last <- subset(d, day==day.last.patient, select=c(taxa.collapsed, count.sample))
colnames(last)[2] <- "label.last"

d <- merge(d, merge(first, last, all=T), all.x=T)

d <- d[!duplicated(d$taxa.collapsed),]
rownames(d) <- d$taxa.collapsed

main <- sapply(1:27, function (patient) {
    # Read the tree, and root at the most prevalent baseline sequence.
    t <- unroot(read.tree(sprintf("06_tree/%02d.nwk", patient)))
    t <- ladderize(root(t, roots[patient]))
    
    # For each tip and internal node, which time points does it appear at?
    tips <- t$tip.label
    tip.data <- data.frame(first=d[tips,"day.min.taxa"] == d[tips,"day.first.patient"],
                           last=d[tips,"day.max.taxa"] == d[tips,"day.last.patient"],
                           tropism=d[tips, "tropism"])
    node.data <- data.frame(t(sapply(prop.part(t), function (clade) {
        tips <- t$tip.label[clade]
        clade.data <- d[tips,]
        r5.count <- sum(clade.data[clade.data$tropism == "R5","count"])
        x4.count <- sum(clade.data[clade.data$tropism == "X4","count"])
        c(first=any(clade.data$day.min.taxa == clade.data$day.first.patient),
           last=any(clade.data$day.max.taxa == clade.data$day.last.patient),
           tropism=if (r5.count < x4.count) "X4" else "R5")
    })))
    tree.data <- rbind(tip.data, node.data)
    
    # Reconstruct tropism at internal nodes.
    #tip.tropism <- d[tips, "tropism"]
    #node.tropism <- MPR(tip.tropism, unroot(t), roots[patient])[,2]
    #tree.data$tropism <- levels(tip.tropism)[c(tip.tropism, node.tropism)]
    
    # Color internal branches by tropism
    tree.data$color <- ifelse(tree.data$tropism == "X4", "red", "blue")
    
    pdf(sprintf("07_pattrees/%02d.pdf", patient), family="Times")
    . <- sapply(c("first", "last"), function (timepoint) {
        t$tip.label <- d[tips, paste0("label.", timepoint)]
        tip.color <- ifelse(tree.data[,timepoint], tree.data$color, "lightgrey")
        edge.color <- tip.color[t$edge[,2]]
        plot(t, type="fan", tip.color=tip.color, edge.color=edge.color)
        mtext(paste("Patient", patient, "day", d[tips[1], paste("day", timepoint, "patient", sep=".")]), side=1, font=2)
    })
    dev.off()
})
