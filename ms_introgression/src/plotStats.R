a <- read.table("stats.txt")

dim(a)

pdf("stats.pdf")
for( i in 1:ncol(a)){
    plot(a[,i])
}
dev.off()
