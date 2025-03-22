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
## looks OK.

## This will also get the scores for each individual position
system.time(
    lp.seq.union <- window.kmer.dist( as.character(lp.seq[[1]]), dinucleotides, 200, ret.flag=1 )
)
##  user  system elapsed 
## 1.790   1.744   3.534 
dim(lp.seq.union$dist)
## [1] 201  16

system.time(
    lp.seq.union <- window.kmer.dist( as.character(lp.seq[[1]]), dinucleotides, 200, ret.flag=0 )
)
##  user  system elapsed 
## 1.198   0.066   1.265 

## a random sequence. This is unfortunately slow; sampling 5 * 48e6 times and pasting.. 
lp.rnd <- paste0( sample(kmer.seq(1), 5 * width(lp.seq)[1], replace=TRUE, prob=lp.mnf$f), collapse="" )
lp.rnd.wc <- window.kmer.dist( lp.rnd, dinucleotides, 200 )
lp.wc <- lp.seq.union

lp.rnd.2 <- paste0( sample(kmer.seq(1), width(lp.seq)[1], replace=TRUE, prob=lp.mnf$f), collapse="" )
lp.rnd.2.wc <- window.kmer.dist( lp.rnd, dinucleotides, 200 )

## do also for monomers:
lp.rnd.mnw <- window.kmer.dist( lp.rnd, kmer.seq(1), 200 )
lp.mnw <- window.kmer.dist( as.character(lp.seq[[1]]), kmer.seq(1), 200 )

dim(lp.wc$dist)
dim(lp.rnd.wc$dist)

lp.rnd.k1 <- kmer.counts( lp.rnd, k=1, with.f=TRUE )
lp.rnd.k2 <- kmer.counts( lp.rnd, k=2, with.f=TRUE )
lp.rnd.dn.exp <- with(lp.rnd.k1, f %*% t(f))
lp.rnd.dn.exp <- as.vector(lp.rnd.dn.exp)
names(lp.rnd.dn.exp) <- names(dn.exp)

## dinucleotide binomial expected counts for a window:
## to compensate for the fact that the distribution is not random...
## ws: window size, dn.p: the expected probability of the dinucleotide in a random
##     sequence.
dn.bin <- function(ws, dn.p, homodimer=FALSE){
    ## for a dinucleotide the effective window size is ws - 1
    ws <- ws - 1
    r1 <- ws %/% 2
    r2 <- r1 + ws %% 2
    ## this considers the sampling of even and odd dinucleotides separately;
    r1.b <- dbinom( 0:r1, r1, dn.p )
    r2.b <- dbinom( 0:r2, r2, dn.p )
    p.m <- r1.b %*% t(rev(r2.b))
    ## this means that r2.b --> columns in reverse order
    n <- (r2 - col(p.m) + 1) + (row(p.m)-1)
    ## if we have homodimer, then we can have up to r1 + r2 instances.
    if(homodimer){
        p <- as.vector(tapply(p.m, n, sum)) / sum(p.m)
        n <- unique(sort(n))
        return(cbind(n, p))
    }
    ## we are still over counting; because selecting a single dinucleotide
    ## at position i, will block positions i-1 and i+1, unless it is adjacent
    ## to a previously selected one.
    ## we can determine the expected number of sites that are next to each other
    ## on the basis of the expected distance. But it gets complicated to do it properly.
    n <- n[ upper.tri(n, diag=FALSE) ]
    p <- p.m[ upper.tri( p.m, diag=FALSE ) ]
    p <- as.vector( tapply( p, n, sum ))
    cbind(unique(sort(n)), p)
}

## ws: window size
## mnf: mononucleotide frequencies
model.window <- function(ws, mnf, nucs){
    dnf <- mnf %*% t(mnf)
    dn <- cbind( nucs[row(dnf)], nucs[col(dnf)] )
    homod <- dn[,1] == dn[,2]
    m <- matrix(nrow=length(dnf), ncol=ws-1)
    trans.m <- sapply(1:16, function(i){ dn[i,2] == dn[,1] }) + 0
    mnf.t <- as.vector( matrix(mnf, nrow=4, ncol=4, byrow=TRUE) )
    m[,1] <- dnf
##    trans.m <- sapply(0:15 %% 4, function(i){ m[,1] * (0:15 %/% 4 == i) })
    for(i in 2:ncol(m)){
        m[,i] <- (trans.m %*% m[,1]) * mnf.t
    }
    m
}

?dpois

ws <- nrow(lp.wc$dist) - 1
par(mfrow=c(4,4))
for(dn in colnames(lp.wc$dist)){
    homodimer <- length( table( unlist(strsplit(dn, ""))) ) == 1
##    n <- ifelse(homodimer, ws-1, ceiling((ws-1)/2))
    n <- ws-1
    bin.1 <- dbinom( 0:n, n, p=dn.exp[dn] )
    bin.2 <- dn.bin( ws-1, dn.exp[dn], homodimer )
    bin.3 <- dbinom( 0:n, n, p=lp.rnd.dn.exp[dn] )
    pois <- dpois( 0:n, lambda=(n * lp.rnd.dn.exp[dn] ))
    plot(1:nrow(lp.wc$dist)-1, lp.wc$dist[,dn], main=dn, type='l', ylim=range(c(lp.wc$dist[,dn], lp.rnd.wc$dist[,dn], bin.1)), lwd=2 )
    lines(1:nrow(lp.rnd.wc$dist)-1, lp.rnd.wc$dist[,dn], col='red', lwd=2)
    lines(1:nrow(lp.rnd.2.wc$dist)-1, lp.rnd.2.wc$dist[,dn], col='red', lwd=2)
    lines(0:n, bin.1, col='blue', lwd=2)
    lines(bin.2[,1], bin.2[,2], col="green", lwd=2)
##    lines(1:nrow(lp.rnd.wc$dist)-1, bin.3, col="purple", lwd=2)
    lines(0:n, pois, col="purple", lwd=2)
}

par(mfrow=c(2,2))
for(i in 1:ncol(lp.mnw$dist)){
    n <- nrow(lp.mnw$dist) - 1
    bin.1 <- dbinom( 0:n, n, p=lp.mnf$f[i] )
    nuc <- colnames(lp.mnw$dist)[i]
    with(lp.mnw, plot(0:n, dist[,i], main=nuc, type='l', ylim=range(c(dist[,i], lp.rnd.mnw$dist[,i], bin.1)), lwd=2, col='black'))
    lines( 0:n, bin.1, col='blue', lwd=2 )
    with(lp.rnd.mnw, lines(0:n, dist[,i], type='l', col='red', lwd=2))
}
### The mononucleotides fit perfectly; And all show a larger spread. A and T have a small peak at 100, suggesting that
### we do have some TA repeats of 200 bases or more. This is consistent with an equal small peak of 0 C or Gs.
### But this does tell us that the sequence generation is correct; but that the distribution methods chosen here are
### not sufficiently good. 

## Lets consider the expected distances between adjacent positions:
gc.pos <- matchPattern( DNAString("GC"), lp.seq[[1]] )
gc.rnd.pos <- start(matchPattern( DNAString("GC"), DNAString(lp.rnd) ))
gc.rnd.pos.e <- gc.rnd.pos[ gc.rnd.pos %% 2 == 0 ]
gc.rnd.pos.o <- gc.rnd.pos[ gc.rnd.pos %% 2 == 1 ]

plot(0:100, dbinom(0:100, 100, 1 - dn.exp['GC']), type='l')
plot(0:100, dbinom(0:100, 100, dn.exp['GC']), type='l')
## all distances will be even
all.dst <- diff(gc.rnd.pos)
even.dst <- diff(gc.rnd.pos.e) / 2
odd.dst <- diff(gc.rnd.pos.o) / 2

all.dst.tb <- as.numeric(table( c(2:max(all.dst), all.dst))) - 1
even.dst.tb <- as.numeric(table( c(1:max(even.dst), even.dst) )) - 1
odd.dst.tb <- as.numeric(table( c(1:max(odd.dst), odd.dst))) - 1

even.dst.tb.f <- even.dst.tb / sum(even.dst.tb)
odd.dst.tb.f <- odd.dst.tb / sum(odd.dst.tb)

even.dst.h <- hist(even.dst, breaks=seq(0.5, max(even.dst)+0.5))
odd.dst.h <- hist(odd.dst, breaks=seq(0.5, max(odd.dst)+0.5))

## expected distribution is
plot( 1:100, (1-dn.exp['GC'])^(1:100) / sum( (1-dn.exp['GC'])^(1:100) ), type='l')
plot( 1:100, even.dst.h$density[1:100], type='l' )
lines( 1:100, (1-dn.exp['GC'])^(1:100) / (1/log(1-dn.exp['GC'])), type='l', col='purple')
lines(1:100, even.dst.tb.f[1:100], col='red')
lines(1:100, odd.dst.tb.f[1:100], col='blue')

## from the internet, the integration of an exponential decay f(x) = exp(-x)
## between points a and b (where b > a ?) is simply (exp(-a) - exp(-b))
## in my case f(x) = exp( log(p)x )

x <- 1:100
p <- 1 - dn.exp['GC']
lp <- log(p)
plot( x, even.dst.h$density[x], type='l' )
lines( x, exp(lp * (x-1)) - exp(lp * x), col='red')

dn.dst <- lapply(names(dn.exp), function(dn){
    rnd <- start(matchPattern( DNAString(dn), DNAString(lp.rnd) ))
    lp <- start(matchPattern( DNAString(dn), lp.seq[[1]] ))
    starts <- list(rnd.even=rnd[ rnd %% 2 == 0 ],
                   rnd.odd=rnd[ rnd %% 2 == 1 ],
                   lp.even=lp[ lp %% 2 == 0 ],
                   lp.odd=lp[ lp %% 2 == 1 ])
    starts.h <- lapply(starts, function(x){
        d <- diff(x/2)
        hist(d, breaks=seq(0.5, max(d)+0.5, 1), plot=FALSE)
    })
    list(rnd=rnd, lp=lp, starts=starts, starts.h=starts.h)
})
names(dn.dst) <- names(dn.exp)

par(mfrow=c(4,4))
for(dn in names(dn.exp)){
    x <- 1:100
    p <- 1 - dn.exp[dn]
    lp <- log(p)
    with(dn.dst[[dn]]$starts.h, {
        plot( x, rnd.even$density[x], type='l', main=dn, col=1, ylim=range(rnd.even$density, lp.even$density), lwd=2 )
        lines( x, lp.even$density[x], type='l', col=2, lwd=2 )
    })
    lines( x, exp(lp * (x-1)) - exp(lp * x), col=4, lwd=2)
}
    
## use a Markov kind of model to estimate probability of a specific dinucleotide at
## all positions in a window:
## p1 and p2; the frequency of nucleotide 1 and 2
## w; the size of the window
hetero.dn.markov.p <- function(p1, p2, w){
    m <- matrix(nrow=w+1, ncol=3)
    ## The columns are the states of the model;
    ## S: starting point. The nucleotide is undefined
    ## n1: At the first nucleotide of the dimer
    ## n2: At the second nucleotide of the dimer
    colnames(m) <- c("S", "n1", "n2")
    m[1,] <- c(1, 0, 0)
    p1.2 <- p1 + p2
    for(i in 2:nrow(m)){
        k <- i-1
        m[i,'S'] <- m[k,'S'] * (1-p1) + m[k,'n1'] * (1-p1.2) + m[k,'n2'] * (1-p1)
        m[i,'n1'] <- sum( p1 * m[k,] )
        m[i,'n2'] <- m[k,'n1'] * p2
    }
    m
}


## lets try for AC
h.dn <- with(lp.mnf, hetero.dn.markov.p(f[1], f[2], 200))

## let's compare...
with(lp.rnd.wc, plot(1:nrow(dist)-1, dist[,'AC'], type='l', lwd=2))
bin.1 <- with(lp.rnd.wc, dbinom( 1:nrow(dist)-1, nrow(dist)-3, mean(dn.exp[c('AC', 'CA', 'GT', 'TG')])))
bin.2 <- with(lp.rnd.wc, dbinom( 1:nrow(dist)-1, nrow(dist)-3, h.dn[201,'n2']))
bin.3 <- with(lp.rnd.wc, dn.bin( nrow(dist)-3, dn.exp['AC'], homodimer ))
lines(1:length(bin.1) - 1, bin.1, col='red', lwd=2)
lines(1:length(bin.2) - 1, bin.2, col='blue', lwd=2)
lines(bin.3[,'n'], bin.3[,'p'], col='gold', lwd=2)
## and we still do not get the numbers that make sense.. 

## this takes a bit of time. 
lp.all.1 <- kmer.counts( as.character(lp.seq), 1 )
lp.all.1$f
## [1] 0.2998180 0.2004253 0.2993487 0.2004079
## 
lp.mnf$f
## [1] 0.3006933 0.1988830 0.3008665 0.1995572

kmer.counts( lp.rnd, 1 )$counts / lp.mnf$counts
## [1] 1.0002718 0.9998203 0.9998656 0.9999826

#### is window.kmer.dist actually giving me correct counts of dinucleotides?
#### I can't seem to find a distribution that describes the random situation
#### so first double check that the function does what it should do.

## This unit contains:
## CG: 2 (one at begin one at at end
## GC: 2 (internal, not repeating)
## CC: 1 internal
## CA: 1 internal
## AA: 1 internal
## AT: 1 internal
## TG: 1 internal
test.monomer <- "CGCCAATGCG"
## check:
data.frame(dn=kmer.seq(2), count=kmer.counts(test.monomer, 2)$counts)

## Aah, this misses one CG; presumably the last one. Check:
data.frame(dn=kmer.seq(2), count=kmer.counts(paste0(test.monomer, "GG"), 2)$counts)
## That gives us 2 CG; but only one GG. Hence the last nucleotide is missing
## from the count.


test.seq.1 <- paste(rep("CGCCAATGCG", 100)


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

