require("Biostrings")
source("kmer_spans.R")

seq <- paste(rep("ATAGACTAATACCTATACTAGACGTACTAGACCGAT", 10), collapse="")
seq.2 <- paste( sample(c("A", "C", "T", "G"), 5e7, replace=TRUE, prob=c(0.3, 0.2, 0.3, 0.2)), collapse="" )

## test windows
seq <- vector(mode='character')
seq[1] <- paste( sample(c("A", "C", "T", "G"), 100, replace=TRUE), collapse="" )
seq[2] <- paste( sample(c("A", "C", "T", "G"), 100, replace=TRUE), collapse="" )
seq[3] <- paste(paste( rep("AG", 50), collapse = ""), seq[1], seq[2], seq[1], sep="")


dyn.load("src/kmer_spans.so")

tmp <- window.kmer.dist(substring(seq[3], 1, 120), c("AG", "GA", "AT"), 20)


tmp <- window.kmer.dist( seq[3], c("A", "C", "T", "G"), 20)


kmer.count <- function(seq, k){
    tmp <- .Call("kmer_counts", seq, as.integer(k))
    names(tmp) <- c("n", "counts")
    c(tmp, k=k)
}

kmer.seq <- function(k){
    .Call("kmer_seq_r", as.integer(k))
}


kmers <- kmer.seq(2)

system.time(
    tmp <- kmer.count(seq.2, 2)
)
names(tmp$counts)  <- kmers

cbind(kmers, tmp$counts)

system.time(
    tmp <- kmer.count(seq.2, 8)
)
##  user  system elapsed 
##  0.045   0.000   0.045 

system.time(
    tmp <- kmer.count(seq.2, 10)
)
##  user  system elapsed 
## 0.133   0.000   0.134 

system.time(
    tmp <- kmer.count(seq.2, 12)
)
##  user  system elapsed
##  0.752   0.024   0.777 

system.time(
    tmp <- kmer.count(seq.2, 6)
)
##  user  system elapsed
## 0.043   0.000   0.042

ns <- "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
nchar(ns)
seq.3 <- paste( seq, ns, seq, sep="" )

dyn.load("src/kmer_spans.so")
tmp2 <- kmer.count(seq, 2)
tmp3 <- kmer.count(seq.3, 2)
tmp4 <- kmer.count(c(seq, seq.3), 2)
cbind(tmp2$counts, tmp3$counts, tmp4$counts)

all(tmp2$counts * 2 == tmp3$counts) ## TRUE
## N's are removed.


## test regions
seq <- vector(mode='character')
seq[1] <- paste( sample(c("A", "C", "T", "G"), 100, replace=TRUE), collapse="" )
seq[2] <- paste( sample(c("A", "C", "T", "G"), 100, replace=TRUE), collapse="" )
seq[3] <- paste(paste( rep("AG", 50), collapse = ""), seq[1], seq[2], seq[1], sep="")




lc.regions <- function(seq, k, min.w, min.score, thr=0.5){
    tmp <- .Call("kmer_low_comp_regions", seq, as.integer(k),
                 as.integer(min.w), as.double(min.score), thr)
    names(tmp) <- c("n", "counts", "w.rank", "pos", "score")
    tmp$pos <- t(tmp$pos)
    tmp$score <- t(tmp$score)
    tmp
}

regs <- lc.regions( seq, 2, 20, 10, 0.5 )
regs.3 <- lc.regions( seq, 3, 20, 10 )
regs.4 <- lc.regions( seq, 4, 20, 10 )
regs.5 <- lc.regions( seq, 5, 20, 10 )

## and lets try with a proper genome:
lp.seq <- readDNAStringSet( "~/genomes/lophius/hifi_asm/yahs.out_scaffolds_final.fa" )
## these are ordered by length. The longest is 48 Mbp.

## lets get the distributions of the lengths:
## for all dinucleotides
dinucleotides <- kmer.seq(2)

## mononucleotide frequencies
lp.mnf <- kmer.counts( as.character(lp.seq[[1]]), 1 )

dn.exp.m <- with(lp.mnf, f %*% t(f))
colnames(dn.exp.m) <- kmer.seq(1)
rownames(dn.exp.m) <- kmer.seq(1)
dn.exp <- as.vector(dn.exp.m)
names(dn.exp) <- paste0( rownames(dn.exp.m)[row(dn.exp.m)], colnames(dn.exp.m)[col(dn.exp.m)] )

## scaffold 11 gives me strange values even within the first million base pairs.
dn.wind.11 <- window.kmer.dist( as.character(lp.seq[[11]]), dinucleotides, 200, freq=FALSE, ret.flag=1L)

rng <- 1:1e7
with(dn.wind.11, plot(rng, scores[[1]][rng,'CG'], type='l'))
## something weird after 7e5
abline(v=6.85e5, col='red')

with(dn.wind.11, plot(1:nrow(scores[[1]]), scores[[1]][,'CG'], type='l', xlim=c(6.8e5, 6.9e5)))
abline(v=c(685700, 685700+500), col='red')
rng <- c(685700, 685700+1000)
with(dn.wind.11, plot(1:nrow(scores[[1]]), scores[[1]][,'CG'], type='l', xlim=c(6.8e5, 6.9e5)))

rng <- (1e6+20):(1e6+400)
with(dn.wind.11, plot(rng, scores[[1]][rng,'CG'], type='l'))

as.character(substring( lp.seq[[11]], rng[1], rng[2]))

with(dn.wind.11, scores[[1]][ 400 + rng[1]:rng[2], 'CG' ])

as.character(substring( lp.seq[[11]], 400 + rng[1], 500 + rng[1] ))

system.time(
    lp.seq.union <- window.kmer.dist( as.character(lp.seq[[1]]), dinucleotides, 200, ret.flag=1 )
)
##  user  system elapsed 
##  1.673   1.149   2.822 
dim(lp.seq.union)
## [1] 201  16

system.time(
    lp.seq.union <- window.kmer.dist( as.character(lp.seq[[1]]), dinucleotides, 200, ret.flag=0 )
)
##  user  system elapsed 
## 1.193   0.064   1.256 

## a random sequence. This is unfortunately slow; sampling 48e6 times and pasting.. 
lp.rnd <- paste0( sample(kmer.seq(1), width(lp.seq)[1], replace=TRUE, prob=lp.mnf$f), collapse="" )
lp.rnd.wc <- window.kmer.dist( lp.rnd, dinucleotides, 200 )

tmp <- kmer.counts( lp.rnd, 1)
tmp$f
## [1] 0.3007744 0.1988468 0.3008254 0.1995534
lp.mnf$f
## [1] 0.3006933 0.1988830 0.3008665 0.1995572

## this takes a bit of time. 
lp.all.1 <- kmer.counts( as.character(lp.seq), 1 )
lp.all.1$f
## [1] 0.2998180 0.2004253 0.2993487 0.2004079
## 
lp.mnf$f
## [1] 0.3006933 0.1988830 0.3008665 0.1995572

kmer.counts( lp.rnd, 1 )$counts / lp.mnf$counts
## [1] 1.0002718 0.9998203 0.9998656 0.9999826

im.col <- hcl.colors(128, "YlOrRd", rev = TRUE)

image(1:nrow(lp.seq.union$dist), 1:ncol(lp.seq.union$dist), log(lp.seq.union$dist), col=im.col, axes=FALSE)
axis(1)
axis(2, at=1:ncol(lp.seq.union), labels=colnames(lp.seq.union), las=2)


## 
wind.exp <- sapply(dn.exp, function(p){ dbinom( 1:nrow(lp.seq.union$dist)-1, as.integer(nrow(lp.seq.union$dist)), p )})

par(mfrow=c(2,2))
image(1:nrow(lp.seq.union$dist), 1:ncol(lp.seq.union$dist), (lp.seq.union$dist), col=im.col, axes=TRUE)
image(1:nrow(lp.rnd.wc$dist), 1:ncol(lp.rnd.wc$dist), (lp.rnd.wc$dist), col=im.col, axes=TRUE)
image(1:nrow(wind.exp), 1:ncol(wind.exp)-1, (wind.exp), col=im.col)

dn.cols <- hsv( 1:ncol(lp.seq.union$dist) / (1.25 * ncol(lp.seq.union$dist)), 1, c(0.5, 0.8) )

par(mfrow=c(2,2))
plot(1:nrow(lp.seq.union$dist)-1, lp.seq.union$dist[,1], ylim=range(log2(1+lp.seq.union$dist)), type='n')
names(dn.cols) <- colnames(lp.seq.union$dist)
for(dn in colnames(lp.seq.union$dist)){
    lines(1:nrow(lp.seq.union$dist) - 1, log2(1+lp.seq.union$dist[,dn]), type='l', col=dn.cols[dn], lwd=2)
##    inpt <- readline("next: ")
}
legend('topright', legend=names(dn.cols), dn.cols, lty=1, lwd=2, col=dn.cols)


par(mfrow=c(4,4))
for(dn in colnames(lp.seq.union$dist)){
    ylim <- range( c(lp.seq.union$dist[,dn], lp.rnd.wc$dist[,dn], wind.exp[,dn]))
    plot(1:nrow(lp.seq.union$dist) - 1, lp.seq.union$dist[,dn] / sum(lp.seq.union$dist[,dn]), type='l', col=1, lwd=2, main=dn, ylim=ylim)
    lines(1:nrow(lp.rnd.wc$dist) - 1, lp.rnd.wc$dist[,dn] / sum(lp.rnd.wc$dist[,dn]), type='l', col=2, lty=1, lwd=2)
    lines(1:nrow(wind.exp) - 1, wind.exp[,dn], type='l', col=4, lty=1, lwd=2)
##    inpt <- readline("next: ")
}

## lets look at all the scaffolds in lophius piscatorius seperately:
wind.s <- 200
sc.wind <- lapply( lp.seq, function(x){
    cat(".")
    seq <- as.character(x)
    cat(":")
    l <- nchar(seq)
    mf <- kmer.counts( seq, 1, with.f=TRUE )
    cat("!")
    wind.f <-  window.kmer.dist( seq, dinucleotides, wind.s, freq=TRUE )
    cat(" ")
    list(l=l, mf=mf, wind=wind.f)
})
length(sc.wind) ## 154

plot( 1:length(sc.wind), sapply(sc.wind, function(x){ x$l }) )
lp.sc.size <- sapply(sc.wind, function(x){ x$l })

min.size <- 1e7
b1 <- lp.sc.size >= min.size
chr.col <- hsv( 0:sum(b1) / (1.2 * sum(b1)), 1, 0.8, 0.4 )
par(mfrow=c(1,1))
for(dn in dinucleotides){
    wf <- sapply(sc.wind[b1], function(x){ x$wind$dist[,dn] })
    col.max <- apply(wf, 2, max)
    plot(1:nrow(wf) - 1, wf[,1], type='n', main=dn, ylim=range(c(wf, lp.rnd.wc$dist[,dn])), xlim=c(0,75))
    x0 <- 2:nrow(wf) - 2
    x1 <- x0 + 1
    y0 <- wf[ 2:nrow(wf)-1, ]
    y1 <- wf[ 2:nrow(wf), ]
    cols <- matrix(chr.col[ col(y0) ], nrow=nrow(y0))
    segments( x0, y0, x1, y1, col=cols, lwd=3)
    lines( 1:nrow(lp.rnd.wc$dist)-1, lp.rnd.wc$dist[,dn], lwd=3, lty=3, col=1)
    legend('topright', legend=paste(names(col.max), sprintf("%.3f", col.max)), lty=1, col=cols[1,], lwd=3 )
    inpt <- readline("next: ")
}
### There are some small differences between the scaffolds:
### in particular, scaffold 3 looks rather different.

tmp <- sapply(sc.wind[b1], function(x){ as.numeric(x$wind$dist) })
sc.pca <- prcomp(t(tmp), scale=FALSE, center=FALSE)

plot(sc.pca) ## basically there's only one dimension here
with(sc.pca, plot(x[,1], x[,2], type='n'))
with(sc.pca, text(x[,1], x[,2], 1:nrow(x)))
## 3 looks different in component 2, not 1. But the range of component 2 is actually
## much greater.

sc.wind.max <- sapply( sc.wind, function(x){
    apply(x$wind$dist, 2, max)
})

sc.wind.max.i <- sapply( sc.wind, function(x){
    apply(x$wind$dist, 2, which.max)
})


par(mfrow=c(1,2))
par(mar=c(5.1, 12.1, 4.1, 2.1))
o <- order( colSums(sc.wind.max)[b1], decreasing=TRUE )
image(x=1:nrow(sc.wind.max), y=1:sum(b1), z=(sc.wind.max[,b1])[,o], col=im.col, axes=FALSE, xlab="", ylab="")
axis(1, at=1:nrow(sc.wind.max), labels=rownames(sc.wind.max))
axis(2, at=1:sum(b1), (colnames(sc.wind.max)[b1])[o], las=2)
##
## par(mfrow=c(1,1))
par(mar=c(5.1, 12.1, 4.1, 2.1))
image(x=1:nrow(sc.wind.max.i), y=1:sum(b1), z=(sqrt(sc.wind.max.i[,b1])[,o]), col=im.col, axes=FALSE, xlab="", ylab="")
axis(1, at=1:nrow(sc.wind.max), labels=rownames(sc.wind.max))
axis(2, at=1:sum(b1), (colnames(sc.wind.max)[b1])[o], las=2)


## will it play nicely with mclapply ?
require(parallel)

wind.s <- 200
system.time(
sc.wind <- mclapply( lp.seq, function(x){
    cat(".")
    seq <- as.character(x)
    cat(":")
    l <- nchar(seq)
    mf <- kmer.counts( seq, 1, with.f=TRUE )
    cat("!")
    wind.f <-  window.kmer.dist( seq, dinucleotides, wind.s, freq=TRUE )
    cat(" ")
    list(l=l, mf=mf, wind=wind.f)
}, mc.cores=20)
)
## system elapsed 
## 23.409   4.619   3.048 

length(sc.wind) ## 154


GCF.seq <- readDNAStringSet( "~/tine/Teleostei_assemblies/teleostei_prot_seq_ncbi/ncbi_dataset/data/GCF_902827115.1/cds_from_genomic.fna" )

sum(nchar(test.seq)) ## 76062661
test.cnt1 <- kmer.count( as.character(test.seq), 2 )
names(test.cnt1$counts) <- kmer.seq(2)
test.cnt2 <- oligonucleotideFrequency( test.seq, width=2 )

test.cnt2.s <- colSums(test.cnt2)

o <- order(kmer.seq(2))
plot( test.cnt1$counts[o], test.cnt2.s )

comp <- sapply( strsplit(names(test.cnt2.s), ""), function(x){
    cmp <- c(A="T", C="G", T="A", G="C")
    paste( cmp[ x ], collapse="" )
})
names(comp) <- names(test.cnt2.s)

plot(test.cnt1$counts, test.cnt1$counts[ comp[names(test.cnt1$counts)] ])

dyn.load("src/kmer_spans.so")
lp.regs.1 <- lc.regions( as.character(lp.seq[1:2]), 6, 100, 20, thr=0.5 )
lp.regs.2 <- lc.regions( as.character(lp.seq[1:2]), 6, 100, 20, thr=0.55 )

par(mfrow=c(1,2))
with(lp.regs.1, hist(log2(pos[,3]-pos[,2])))
with(lp.regs.2, hist(log2(pos[,3]-pos[,2])))

with(lp.regs.1, plot(sort(log2(pos[,3]-pos[,2]))))
with(lp.regs.2, plot(sort(log2(pos[,3]-pos[,2]))))

with(lp.regs.1, sum(pos[,3] - pos[,2] > 1000) / nrow(pos) ) ## 0.2752913
with(lp.regs.2, sum(pos[,3] - pos[,2] > 1000) / nrow(pos) ) ## 0.1361541

with(lp.regs.1, sum(pos[,3] - pos[,2]))/1e6 ## 80 million?
with(lp.regs.2, sum(pos[,3] - pos[,2]))/1e6 ## 19 million?

plot.regions <- function(i, regs, tform=log10){
    plot.new()
    b <- regs$pos[,1] == i
    scores <- tform(regs$score[b,1])
    plot.window(xlim=c(0, nchar(lp.seq[i+1])), ylim=range(scores))
    with(regs, segments(pos[b,2], scores, pos[b,3], scores))
    axis(1)
    axis(2)
}

par(mfrow=c(2,1))
plot.regions(0, lp.regs.1)
plot.regions(0, lp.regs.2)

par(mfrow=c(2,1))
plot.regions(1, lp.regs.1)
plot.regions(1, lp.regs.2)


lp.regs.3 <- lc.regions( as.character(lp.seq[1:2]), 10, 100, 20, thr=0.5 )
lp.regs.4 <- lc.regions( as.character(lp.seq[1:2]), 10, 100, 20, thr=0.55 )

lp.regs.4.2 <- lc.regions( as.character(lp.seq[1:23]), 10, 100, 50, thr=0.6 )

par(mfrow=c(2,1))
plot.regions(0, lp.regs.3)
plot.regions(0, lp.regs.4)

par(mfrow=c(2,1))
plot.regions(1, lp.regs.3)
plot.regions(1, lp.regs.4)

lp.regs.5 <- lc.regions( as.character(lp.seq[1:2]), 12, 100, 20, thr=0.5 )
lp.regs.6 <- lc.regions( as.character(lp.seq[1:2]), 12, 100, 20, thr=0.55 )

par(mfrow=c(2,1))
plot.regions(0, lp.regs.5)
plot.regions(0, lp.regs.6)

par(mfrow=c(2,1))
plot.regions(1, lp.regs.5)
plot.regions(1, lp.regs.6)

par(mfrow=c(3, 2))
plot.regions(0, lp.regs.1)
plot.regions(0, lp.regs.2)
plot.regions(0, lp.regs.3)
plot.regions(0, lp.regs.4)
plot.regions(0, lp.regs.5)
plot.regions(0, lp.regs.6)

par(mfrow=c(3, 2))
plot.regions(1, lp.regs.1)
plot.regions(1, lp.regs.2)
plot.regions(1, lp.regs.3)
plot.regions(1, lp.regs.4)
plot.regions(1, lp.regs.5)
plot.regions(1, lp.regs.6)


k <- c(8, 10, 12, 14)
## there are essentially 23 chromosomes
lp.regs <- lapply(k, function(x){
    lc.regions( as.character(lp.seq[1:23]), x, 100, 50, thr=0.6 )
})
               

with(lp.regs.2, hist(log2(pos[,3]-pos[,2])))

sapply(lp.regs, function(x){ dim(x$pos) })
##       [,1] [,2]  [,3]  [,4]
## [1,] 24611    0 45268 51972
## [2,]     3    3     3     3
##
## no windows at all with k=10 and thr=0.6
## that's rather strange?

## plot by chromosome the value from each:
par(mfrow=c(6,4))
ylim <- log10(range(unlist(sapply(lp.regs, function(x){ x$score[,1] }))))
for(chr in 0:22){
    xlim <- c(1, nchar(lp.seq[chr+1]))
    plot(1,1, type='n', xlim=xlim, ylim=ylim, xlab="pos", ylab="score")
    for(j in 1:length(lp.regs)){
        b <- lp.regs[[j]]$pos[,1] == chr;
        with(lp.regs[[j]], segments( pos[b,2], log10(score[b,1]), pos[b,3], log10(score[b,1]), col=j))
    }
}
    
par(mfrow=c(2,2))
for(i in c(1,3,4)){
    with(lp.regs[[i]], hist( log10( pos[,3]-pos[,2])))
}

## and lets try to extract the sequences and determine kmers;
## count tetramers
reg.kmer.sp <- lapply(c(1,3,4), function(i){
    ss <- with(lp.regs[[i]], subseq( lp.seq[ pos[,1]+1], pos[,2], pos[,3] ))
    fwd <- oligonucleotideFrequency( ss, width=4, as.prob=TRUE )
    rev <- oligonucleotideFrequency( reverseComplement( ss ), width=4, as.prob=TRUE )
    (fwd + rev) / 2
})

reg.seq <- lapply(c(1,3,4), function(i){
    with(lp.regs[[i]], subseq( lp.seq[ pos[,1]+1], pos[,2], pos[,3] ))
})

reg.ent <- lapply(reg.kmer.sp, function(x){
    apply(x, 1, function(y){
        -1 * sum( ifelse(y > 0, log2(y), 0) * y )
    })
})

par(mfrow=c(2,2))
for(x in reg.ent){
    plot(sort(x), type='l')
}

par(mfrow=c(2,2))
for(i in 1:length(reg.seq)){
    plot(log10(nchar(reg.seq[[i]])), reg.ent[[i]], cex=0.5, col=rgb(0,0,0,0.2))
}


## that's very quick.. ;-)
reg.kmer.sp.pca <- lapply(reg.kmer.sp, function(x){
    prcomp( x, center=FALSE )
})
## but this is really slow!! 

par(mfrow=c(2,2))
for(x in reg.kmer.sp.pca)
    plot(x)
## three dimensions needed to define, but what about..

par(mfrow=c(2,2))
for(p in reg.kmer.sp.pca)
    with(p, plot(x[,1], x[,2]))

par(mfrow=c(2,2))
with(reg.kmer.sp.pca[[1]], plot(x[,1], x[,2]))
with(reg.kmer.sp.pca[[1]], plot(x[,1], x[,3]))
with(reg.kmer.sp.pca[[1]], plot(x[,2], x[,3]))

## it seems that most sequences are really close to the center
## we can consider to order by the rowSums of x
o <- order(rowSums(reg.kmer.sp.pca[[1]]$x[,1:3]^2))
tmp <- subseq( lp.seq[ 1 + lp.regs[[1]]$pos[o,1] ], lp.regs[[1]]$pos[o,2], lp.regs[[1]]$pos[o,3])

par(mfrow=c(1,1))
n <- 1000
image(1:256, 1:n, t(reg.kmer.sp[[1]][o[1:n], ]) )

## we can also determine distances and make a tree from those
## ?dist defaults to euclidean distances between rows

## this is slow; reasonably so as it will compute
## 24611^2 distances;
## 
reg.d1 <- dist(reg.kmer.sp[[1]])

## we can make a tree using the hclust function
## or the ape nj function, or upgma? upgma ought to be
## faster, and given that this is pretty slow, it would
## make sense to do something to increase the speed.

## nj, far, far to slow..
reg.hc <- hclust(reg.d1)


plot(reg.hc, labels=FALSE)

## we can convert the hclust object into a phylo object
## which I'm more familiar with.. but first we can try:
par(mfrow=c(1,1))
image( 1:256, 1:length(reg.hc$order), t(reg.kmer.sp[[1]][ reg.hc$order, ]) )

image( 1:length(reg.hc$order), 1:256, reg.kmer.sp[[1]][ reg.hc$order, ] )

plot(1:length(reg.hc$order), reg.ent[[1]][reg.hc$order], type='l')

## have a look at some sequences...
as.character(reg.seq[[1]][ reg.hc$order[10000:10010] ])

## Restrict ourselves to an interesting subset of sequences:
## the longer kmers have much more information:
sapply( reg.ent, function(x){ sum(x > 7) / length(x) })
## [1] 0.1181179 0.4215782 0.5189910

## lets look a the > 7.5 in the last set.. 
sapply( reg.ent, function(x){ sum(x > 7.5) })
## [1]     2  7807 11299

b <- reg.ent[[3]] > 7.7
reg.d2.i <- which(b)
reg.d2 <- dist(reg.kmer.sp[[3]][b, ])
reg.hc2 <- hclust(reg.d2)

plot(reg.hc2, labels=FALSE)


## lets see if a PCA gives us something more interesting:
reg.pca <- prcomp( reg.kmer.sp[[3]][b, ] )
## very little structure in that..
with(reg.pca, plot(x[,1], x[,2]))

o.i <- reg.d2.i[ reg.hc2$order ]
image(1:length(o.i), 1:256, reg.kmer.sp[[3]][o.i, ] )

## we can also do:
i <- reg.hc2$order
image( as.matrix(reg.d2)[i,i] )

require(ape)

reg.ph <- as.phylo( reg.hc2 )

## traverse an ape tree and get groups where the max distance to the
## leaf from the current node is max.dist
get.ph.groups <- function(tree, max.dist){
    descend <- function(root){
        root.i <- which(tree$edge[,1] == root)
##        browser()
        if(length(root.i) == 0){
            i <- which(tree$edge[,2] == root)
            return(matrix(c(i, root, 0),
                          nrow=1))
        }
        ## the length is added on here, in order to be cumulative.
        children <- lapply(root.i, function(r){
            tmp <- descend(tree$edge[r, 2])
            tmp[,3] <- tmp[,3] + tree$edge.length[r]
            tmp
        })
        names(children) <- tree$edge[root.i, 2]
        max.d <- max(sapply(children, function(x){ max(x[,3]) }))
        if(max.d <= max.dist)
            return( do.call(rbind, children) )
        node.groups <<- c(node.groups, children)
        return(matrix(nrow=0, ncol=3))
    }
    root <- setdiff( tree$edge[,1], tree$edge[,2] )
    node.groups <- list()
    descend(root)
    return(node.groups)
}

tmp <- get.ph.groups( reg.ph, 0.015 )

## and this takes for ever. Which I guess is not that strange
## given that a lot of distance matrices will need to be calculated

## this is way too slow.
## system.time(
##     reg.nj <- nj(reg.d1)
## )

with(lp.regs, plot( (pos[,3] - pos[,2]), score[,1] ))
## the log log is highly linear; 
with(lp.regs, plot( log10(pos[,3] - pos[,2]), log10(score[,1]) ))

plot(log10(lp.regs$score[,1]))

hist( log10(lp.regs$score) )
with(lp.regs, hist( log2(pos[,3] - pos[,2]), breaks=50))

## we have some bugs:
i <- which(lp.regs$score[,1] < 20)
head( lp.regs$pos[i-10, ], n=20)
head( lp.regs$score[i-10, ], n=20)

o <- order(lp.regs$pos[,3] - lp.regs$pos[,2], decreasing=TRUE)
head(lp.regs$pos[o,])

with(lp.regs, {b=pos[,1] == 0; plot(pos[b,2], score[b,1])})
with(lp.regs, {b=pos[,1] == 1; plot(pos[b,2], score[b,1])})

