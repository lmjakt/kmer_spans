require("Biostrings")
dyn.load("src/kmer_spans.so")

seq <- paste(rep("ATAGACTAATACCTATACTAGACGTACTAGACCGAT", 10), collapse="")
seq.2 <- paste( sample(c("A", "C", "T", "G"), 5e7, replace=TRUE, prob=c(0.3, 0.2, 0.3, 0.2)), collapse="" )

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
## lp.seq <- readDNAStringSet( "~/genomes/lophius/hifi_asm/yahs.out_scaffolds_final.fa" )
## these are ordered by length. The longest is 48 Mbp.

test.seq <- readDNAStringSet( "~/tine/Teleostei_assemblies/teleostei_prot_seq_ncbi/ncbi_dataset/data/GCF_902827115.1/cds_from_genomic.fna" )

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

