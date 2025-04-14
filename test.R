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
lp.mnf <- kmer.counts( as.character(lp.seq[[1]]), 1, with.names=TRUE )
dn.exp.m <- with(lp.mnf, f %*% t(f))
colnames(dn.exp.m) <- kmer.seq(1)
rownames(dn.exp.m) <- kmer.seq(1)
dn.exp <- as.vector(dn.exp.m)
names(dn.exp) <- paste0( rownames(dn.exp.m)[row(dn.exp.m)], colnames(dn.exp.m)[col(dn.exp.m)] )


## This will also get the scores for each individual position
system.time(
    lp.wc <- window.kmer.dist( as.character(lp.seq[[1]]), dinucleotides, 200, ret.flag=0 )
)
##  user  system elapsed 
## 1.790   1.744   3.534 
dim(lp.wc$dist)
## [1] 201  16

require("Biostrings")
source("kmer_spans.R")

## a random sequence. This is unfortunately slow; sampling 20 * 48e6 times and pasting.. 
lp.rnd <- paste0( sample(kmer.seq(1), 20 * width(lp.seq)[1], replace=TRUE, prob=lp.mnf$f), collapse="" )
lp.rnd.wc <- window.kmer.dist( lp.rnd, dinucleotides, 200 )
lp.rnd.wc.100 <- window.kmer.dist( lp.rnd, dinucleotides, 100 )
lp.rnd.wc.2 <- window.kmer.dist( lp.rnd, dinucleotides, 2 )

lp.wc.disc <- window.kmer.dist(as.character(lp.seq[[1]]), dinucleotides, 200, ret.flag=3L )


system.time(
    lp.rnd.wc.disc <- window.kmer.dist(lp.rnd, dinucleotides, 200, ret.flag=3L )
)
## This is now 979.0459 million bases. 
##   user  system elapsed 
## 74.010   1.864  75.876 
## not super fast. but OK.

## this is with discrete windows
## should not differ very much.
lp.rnd.wc.wdisc <- window.kmer.dist(lp.rnd, dinucleotides, 200, ret.flag=4L)



## see if these fit the expected data:
par(mfrow=c(4,4))
for(dn in colnames(lp.wc.disc$dist)){
    x <- 1:nrow(lp.wc.disc$dist) - 1
    exp.f <- dbinom( 0:max(x), size=max(x), prob=dn.exp[dn] )
    x <- 0:40
    ylim <- range(c(exp.f, lp.wc.disc$dist[x+1,dn]))
    plot(x, lp.rnd.wc.disc$dist[x+1,dn], type='l', lwd=2, col="blue", main=dn, ylim=ylim)
    lines(x, lp.wc.disc$dist[x+1,dn], type='l', lwd=2)
    lines(x, exp.f[x+1], col='red', lwd=2)
##    inpt <- readline("next: ")
}
## these fit perfectly. So the counting does seem to be correct. But it would still be good
## to have an expected distribution for the overlapping counts.

## lets divide by GC content of dinucleotides:
## OK, this is a bit of a bug: it doesn't count the last position... 
gc.n <- sapply(dinucleotides, function(x){
    sum( window.kmer.dist(x, c("G", "C"), window=2, freq=FALSE)$dist[ 2:3, ] * 1:2 )
})

par(mfrow=c(2,2))
invisible(tapply(dinucleotides, gc.n, function(dn){
    i <- 1:nrow(lp.wc.disc$dist)
    exp.f <- sapply(dn.exp[dn], function(p){ dbinom(i-1, size=max(i), p) })
    rnd.f <- lp.rnd.wc.disc$dist[,dn]
    lp.f <- lp.wc.disc$dist[,dn]
    i <- 1:40
    x <- i-1
    plot(1, 1, type='n', xlim=range(x), ylim=range(exp.f, rnd.f, lp.f), main=paste(dn, collapse=", "))
    cols <- hsv(1:length(dn) / length(dn), 0.9, c(0.5, 0.95))
    for(j in 1:length(dn)){
        lines(x, lp.f[i,j], col=cols[j], lty=1, lwd=2)
        lines(x, rnd.f[i,j], col=cols[j], lty=2, lwd=2)
        lines(x, exp.f[i,j], col=cols[j], lty=3, lwd=2)
    }
    legend('topright', dn, col=cols, lty=1, lwd=2)
}))


par(mfrow=c(2,2))
exp.rnd.f <- invisible(tapply(dinucleotides, gc.n, function(dn){
    i <- 1:nrow(lp.wc$dist)
    exp.f <- sapply(dn.exp[dn], function(p){ dbinom(i-1, size=max(i), p) })
    rnd.f <- lp.rnd.wc$dist[,dn]
    lp.f <- lp.wc$dist[,dn]
    if(FALSE){
        exp.f <- log(exp.f)
        rnd.f <- log(rnd.f)
        lp.f <- log(lp.f)
    }
    i <- 1:60
    x <- i-1
    plot(1, 1, type='n', xlim=range(x), ylim=range(cbind(exp.f, rnd.f, lp.f)[i,], finite=TRUE), main=paste(dn, collapse=", "))
    cols <- hsv(1:length(dn) / length(dn), 0.9, c(0.5, 0.95))
    for(j in 1:length(dn)){
        lines(x, lp.f[i,j], col=cols[j], lty=1, lwd=2)
        lines(x, rnd.f[i,j], col=cols[j], lty=2, lwd=2)
        lines(x, exp.f[i,j], col=cols[j], lty=3, lwd=2)
    }
    legend('topright', dn, col=cols, lty=1, lwd=2)
    het.b <- sapply( strsplit(dn, ""), function(x){ length(unique(x)) }) == 2
    data.frame( cbind(rnd.f[,het.b], exp.f[,het.b]))
}))

exp.rnd.f.lr <- lapply(exp.rnd.f, function(x){ b <- grepl("1", colnames(x)); log(x[,!b] / x[,b])})
exp.rnd.dn <- sapply(exp.rnd.f.lr, function(x){ dn.exp[ colnames(x) ]})

exp.rnd.f.lr <- as.data.frame(do.call(cbind, exp.rnd.f.lr))
exp.rnd.dn <- unlist(exp.rnd.dn)
colnames(exp.rnd.f.lr) <- sub("^[0-9].", "", colnames(exp.rnd.f.lr))
names(exp.rnd.dn) <- sub("^[0-9].", "", names(exp.rnd.dn))

## then we can plot and see if we can work out a relationship:
i <- 1:60
par(mfrow=c(1,1))
cols <- hsv(1:ncol(exp.rnd.f.lr)/ncol(exp.rnd.f.lr), 1, c(0.5, 0.9))
ylim <- range(exp.rnd.f.lr, finite=TRUE)
ylim <- c(-2.3, 0.12)
plot(1, 1, type='n', xlab="n", ylab="log ratio", xlim=range(i-1), ylim=ylim)
invisible(mapply(function(y, col){ lines(i-1, y[i], col=col, lwd=2) }, exp.rnd.f.lr, cols ))
abline(v=apply( exp.rnd.f.lr, 2, which.max ) - 1, col='red')
abline(v=exp.rnd.dn * 200, col='blue')

## can I fit a polynomial model to this? I suspect that it might be a bit difficult.
## but I can try:
exp.rnd.f.lr.mod <- apply(exp.rnd.f.lr, 2, function(x){
    ## exclude infinite values;
    b <- is.finite(x) & !is.na(x)
    beg <- which(b)[2]
    end <- which(!b[beg:length(x)])[1] - 1
    n <- 1:length(x) - 1
    r1 <- n[beg:end]
    r2 <- r1^2
    r3 <- r1^3
    r4 <- r1^4
    lr <- log(r1)
    lm( x[beg:end] ~ r1 + r2 + r3 + r4 + lr)
##    lm( x[beg:end] ~ r1 + lr)
})

## and that does seem to give nicely fitting values.
tmp <- sapply(exp.rnd.f.lr.mod, function(x){ summary(x)$coefficients[,1] })
par(mfrow=c(2,3))
for(i in 1:nrow(tmp)){
    plot(dn.exp[ colnames(tmp) ], tmp[i,], main=rownames(tmp)[i])
}

i <- 1:60
par(mfrow=c(1,1))
cols <- hsv(1:ncol(exp.rnd.f.lr)/ncol(exp.rnd.f.lr), 1, c(0.5, 0.9))
ylim <- range(exp.rnd.f.lr, finite=TRUE)
ylim <- c(-2.3, 0.12)
plot(1, 1, type='n', xlab="n", ylab="log ratio", xlim=range(i-1), ylim=ylim)
invisible(mapply(function(y, col){ lines(i-1, y[i], col=col, lwd=2) }, exp.rnd.f.lr, cols ))
abline(v=apply( exp.rnd.f.lr, 2, which.max ) - 1, col='red')
abline(v=exp.rnd.dn * 200, col='blue')

for(i in 1:ncol(tmp)){
    x <- 1:40
    cr <- tmp[,i]
    y <- sapply(x, function(z){ cr[1] + sum(cr[-1] * c(z, z^2, z^3, z^4, log(z))) })
    lines( x, y, lwd=4, col="gold" )
    inpt <- readline(paste("next", i))
}


## try to see if we can estimate an effective sampling size for each number
## that is try to find an m, and an M (where mr <- M-m)
## such that f^n * f^mr * choose(M, m) = observed frequency...
mmx <- lapply(1:length(exp.rnd.dn), function(i){
    dn <- names(exp.rnd.dn)[i]
    j <- 1:40
    f.obs <- lp.rnd.wc$dist[j,dn]
    n <- j - 1
##    M <- seq(196, 215, 0.05)
    M <- 150:210
    p <- exp.rnd.dn[i]
    ## it should be possible to vectorise this and express it as
    ## matrix multiplications, but keep it simple..
    t(rbind(f.obs, sapply(1:length(f.obs), function(k){
        ##        dbinom(n[k], M, p)
##        p^(n[k]) * (1-p)^(M-n[k]) * choose(200, n[k])
        p^(n[k]) * (1-p)^(200-n[k]) * choose(M, n[k])
    })))
})

par(mfrow=c(3,4))
n <- ncol(mmx[[i]])
cols <- hsv(1:n/n, 1, (n+1:n)/(2.1 * n))
for(i in 1:length(mmx)){
    j <- 1:nrow(mmx[[i]])
    x <- j - 1
    dn <- names(exp.rnd.dn)[i]
    plot(x, log(lp.rnd.wc$dist[j,dn]), type='l', lwd=3, ylim=range(log(mmx[[i]]), finite=TRUE), main=dn)
    for(k in 1:n){
        lines(x, log(mmx[[i]][j,k]), col=cols[k])
    }
    lines(x, log(lp.rnd.wc$dist[j,dn]), type='l', lwd=3, col='gold')
}

mmx.o <- lapply(mmx, function(mm){
    apply(mm, 1, function(x){
        order(abs(log(x[1] / x[-1])))
    })
})

im.col <- hcl.colors(128, "YlOrRd", rev = TRUE)

par(mfrow=c(3,4))
for(i in 1:length(mmx.o)){
##    image(x=150:210, y=0:39, z=mmx.o[[i]], col=im.col)
    plot(1:ncol(mmx.o[[i]])-1, mmx.o[[i]][1,], type='l')
##    image(x=seq(196, 210, 0.05 ), y=0:39, z=mmx.o[[i]], col=im.col)
}
## this doesn't really help that much. Just getting the correction factor from the
## earlier plots seems more reasonable.. 


mmx.fit <- sapply(1:length(mmx.o), function(i){
    p <- exp.rnd.dn[i]
    n <- 1:ncol(mmx.o[[i]]) - 1
    M <- 100 + mmx.o[[i]][1,] - 1
    p^n * (1-p)^(M-n) * choose(M, n)
})

colnames(mmx.fit) <- names(exp.rnd.dn)
par(mfrow=c(3,4))
for(dn in colnames(mmx.fit)){
    i <- 1:nrow(mmx.fit)
    x <- i-1
    plot(x, lp.rnd.wc$dist[i,dn], type='l', col='black')
    lines(x, mmx.fit[,dn], col='red')
}

##lp.rnd.2 <- paste0( sample(kmer.seq(1), width(lp.seq)[1], replace=TRUE, prob=lp.mnf$f), collapse="" )
lp.rnd.2 <- paste0( sample(kmer.seq(1), 1000, replace=TRUE, prob=lp.mnf$f), collapse="" )
lp.rnd.2.wc <- window.kmer.dist( lp.rnd, dinucleotides, 200, freq=FALSE )

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

## This attempts to account for the chain dependency of k-mer
## identities. 
## The probability of a dimer at position i (p[dm,i]):
## can be estimated as: p[n1] * p[n2]
## however, if the dimer at position i-1 was dm, then p[dm,i] is 0
## and this needs to be taken into account. The first dinucleotide (dm)
## is independent so we can state:
## p[dm,1] = p[n1] * p[n2]
## Then we can define p[dm,i] for i > 2 as:
## p[dm,i] = (p[n1] * p[n2]) * (1 - p[dm,i-1])
##
## That implies that the probability of hetero-dimers decreases the longer
## the window. That is non-intuitive, but we'll see if it works.
## For homo-dimers, the probability is a bit different. If the last one
## was a dimer, then, we have a higher probability. I don't think that
## this should change anything, but we can define similarly:
##
## p[dm,i] = (p[n2] * p[dm,i-1]) + (1 - p[dm,i-1]) * (1-p[n1]) * p[n2]
## 
## If the previous one was dm, then the probability is just p[n2]
## Otherwise, it's the probability of the previous dinucleotide being
## compatible without being dm times the probablity of p[n2].
## at the second dimer position this will be:
## p[n2] * p[n1] * p[n2] + (1 - p[n1] * p[n2]) * (1-p[n1]) * p[n2]
## 
## but here p[n1] == p[n2], so it is simply
## p[n1] * p[n1] * p[n1] + (1 - p[n1] * p[n1]) * (1-p[n1]) * p[n1]
## p[n1] * p[n1] * p[n1] + (1 - p[n1] * p[n1]) * (p[n1]-p[n1]^2)
## p[n1]^3 + (1 - p[n1]^2) * (p[n1]-p[n1]^2) 
## p[n1]^3 + p[n1] - p[n1]^2 -p[n1]^3 + p[n1]^4
## p[n1] - p[n1]^2 + p[n1]^4
## ## but that doesn't seem reasonable, so I must have made a mistake
## of some sort
## 
## ws: the window size
## mnf: the mononucleotide frequencies (named vector)
## nucs: the nucleotides in the dimer
##       in theory, the idea used here could be used for
##       any k-mer; but first check if the idea actually
##       works for dimers
## The function will return a set of dinucleotide probabilities
## that can be used with the Poisson Binomial distribution
model.window.2 <- function(ws, mnf, nucs){
    f1 <- mnf[nucs[1]]
    f2 <- mnf[nucs[2]]
    dn.f <- mnf %*% t(mnf)
    rownames(dn.f) <- names(mnf)
    homodimer <- (nucs[1] == nucs[2])
    p1 <- f1 * f2
    p <- rep(0, ws)
    dm1.compat <- dn.f[ , nucs[2] ]
    ## leave probs[1] as 0, as it represents the first position
    ## in the window; The 2 below could be k
    p[2] <- p1
    for(i in 3:ws){
        if(homodimer)
            p[i] <- (p[i-1] * f2) + (1-p[i-1]) * (1-f1) * f2 * f2
        else
            p[i] <- (1 - p[i-1]) * (f1/(1-p1)) * f2
        ## I believe this is correct; but that actually gives us exactly the same
        ## probs; as that reduces to f1 * f2. But we know that doesn't give use
        ## distributions that fit the random observed one.
##            p[i] <- (1 - p[i-1]) * (f1/(1-f2)) * f2
##            p[i] <- (1 - p[i-1]) * (sum(dm1.compat) / (1-p1)) * f2
##            p[i] <- p1 * (1 - p[i-1])
    }
    p
}

## This considers how the effective number of positions that
## can be selected decreases with the number of positive
## dimers. It is based on the fact for for heterodimers,
## is positive, the two overlapping dimers can no longer
## be selected.
## For heterodimers, we consider that the probability of
## i positions in the window being the dimer is:
##
## p[i] <- (dn.f^i + (1-dn.f)^i) * choose(w, i)
##
## For heterodimers, the maximum number of matching positions
## is (w / 2), since we cannot have a match at position i, and i+1
## However, at the beginning of the selection we do have w positions
## to choose from. The probability of choosing exactly w/2
## is clearly f^(w/2), but the probability of choosing one is
## f + (1-f)^w * choose(w, 1). Hence The effect window size decreases
## with the number of matching positions. To compensate for thiw
## we can define the effective window size (w.e) as:
##
## w.e[i] <- w - (i) * (r[i] + 2*(1-r[i]))
##
## where r[i] gives the expected number of sites that are not singletons
## (distance to closest neighbour is 1)
## and thus block neighbours on side only. For this we consider that
## the window borders are also sites since a match at the beginning or
##
## r[i] <- 2 * (1-f)
## where f is (i+2)/(w + 2).
## (1-f) is the expected number of matched sites where the site
## p[i] - p[i-1] are expected to be 1. We double this since the
## neighbour can be both upstream and downstream.
### Scratch that, for dimers it doesn't matter; r[i] is just 2 * i
### Hence
### w.e[i] <- w.e - 2 * i
### but I may need to define this recursively.
## We can then estimate the likelihood of a given number of matches
model.window.3 <- function(ws, dn.p, homodimer=FALSE){
    ## forward (f.l) and reverse (r.l) likelilhoods
    f.l <- dn.p
    r.ll <- log(1-dn.p)
    ## n is the number of matched sites
    n <- 0:ws
    ## p is the probability of selecting n
    p <- rep(0, 1 + ws)
    ## we need to determine the number of selectable sites.
    ## This depends on the expected distances between sites
    ## where this is given as a
    ##
    ## d.exp <- exp(lp * (d-1)) - exp(lp * d)
    ## d is the distance between sites, and p
    ## is log(p) and p is (n+1)/(w+1).
    ## The addition of 1 is because we consider the 0th index to be a site
    ## the distance expect matrix is:
    ## from the internet, the integration of an exponential decay f(x) = exp(-x)
    ## between points a and b (where b > a ?) is simply (exp(-a) - exp(-b))
    ## in my case f(x) = exp( log(p)x )
    d.e <- sapply( 0:ws, function(i){
        lp <- log(1 - (i+1) / (ws+1) )
        d <- 1:ws
        exp(lp * (d-1)) - exp(lp * d)
    })
    ## d.e[1,] is for adjacency. seperation = 0
    ## d.e[2,] is for a separation of 1.
    ## Which is not possible; we need to remove:
    ## The number of sites that are reduced depends on the distances:
    ## adjacent: 2 sites removed
    ## separated by an even distance: 3 sites   |__**__**__|
    ## separated by an odd distance: 4 sites    |_**___**__|
    ## but that doesn't matter until they fill up. So we will ignore.
    w.e <- sapply(1:ws, function(i){
        if(i < 2)
            return(ws);
        sep.p <- d.e[,i+1]
##        sep.p[2] <- 0
##        sep.p <- sep.p/sum(sep.p)
        ws - (i-1) * (2 * sep.p[1] + (3 * (1 - sep.p[1])))
    })
    ## define the probalities of picking the 1st one, then the second one and
    ## so on. This does not give the probability of 0, but that is easy to get
    p.inc <- 1 - dbinom( 0, as.integer(w.e), dn.p )
    p.inc
    ## ## considers a window of w + k-1, allowing a count of w k-mers.
    ## w.e <- ws - c(0, (1:ws-1) * 2)
    ## j <- which(w.e > 0)
    ## ## p[b] gives the probability of selecting exactly one from the
    ## ## reduced set. There are exactly w.e ways of selecting one
    ## p[j] <- f.l * exp(r.ll * (w.e[j] - 1)) * w.e[j]
    ## p
    ## for(i in 1:ws){
    ##     if(w.e[i] < i)
    ##         break
    ##     p[i] <- (exp(f.l * i) + exp(r.l * (w.e[i] - i))) * choose(w.e[i], i)
    ## }
    ## i <- 0:ws
    ## ## f is the frequency of matching sites
    ## f <- (i+2) / (ws+2)
    ## r <- 2 * (1-f)
    ## ## we will refer the the effective size as w
    ## w <- ws - i * (r + 2*(1-r))
    ## ## consider log transforming for range:
    ## fp <- log(dn.p) ## forward prob
    ## rp <- log(1-dn.p) ## reverse prob
    ## dn.rnd <- rep(0, ws + 1)
    ## dn.rnd <- exp(fp * i) + exp(rp * ws-i) * choose( as.integer(w), i )
}

### model 4
### determine likelihoods of gap lengths, and for each one
### recurse down to the next level.
model.4 <- function(ws, dn.p, homodimer){
    gpl <- log(1-dn.p)
    gap.p <- function(dst){
        exp( gpl * (dst-1)) - exp(gpl * dst)
    }
    lapply(1:length(ws), function(i){
        dst <- (i+1):ws - i
        dst.p <- gap.p(dst)
        ## each position in dst.p gives the likelihood 
        ## a dn in any position there or after. 
        lapply(1:length(dst.p), function(j){
            

}

ws <- 200
at.dp <- model.window.2( ws, lp.mnf$f, c("A", "T") )
cg.dp <- model.window.2( ws, lp.mnf$f, c("C", "G") )
cc.dp <- model.window.2( ws, lp.mnf$f, c("C", "C") )

require(poisbinom)

plot(dpoisbinom(0:ws, cg.dp), type='l')
lines(dbinom(0:ws, ws, dn.exp['CG']), col='red')
lines(dpoisbinom(0:ws, cc.dp), type='l', col='blue')
lines(dbinom(0:ws, ws, dn.exp['CC']), col='purple')
## 

mod5 <- function(p, ws){
    ## if we select m sites from n; how many of those
    ## are illegal selections? It seems that the number
    ## of ways of choosing  should be reduced by this number.
    ## k1, k2, refer to the first kmer found in the set
    ##  m
    ##  1  ws positions; no interference
    ##  2  k1: k1 has ws positions, k2 has ws-1 of which ws-3 interfere
    ##  3  k3 has ws-2 possible positions, of which ws-4 interfere.
    ##  m  m4 has ws-[m-1] possible positionf of which ws-m interfere.
    intf.n <- sapply( 0:ws, function(m){
        if(m < 2)
            return(0)
        sum(sapply(2:m, function(i){
            ws - (i+1)
        }))
    })
    sapply(0:ws, function(i){
        nn <- ws - (i-1) * 2
        nn2 <- ws - (i-1) * 1
        nn <- ifelse(nn > ws, ws, nn)
        c(i, nn,
          exp(log(1-p) * (nn-i) + log(p) * i) * choose(nn, i),
          exp(log(1-p) * (nn2-i) + log(p) * i) * choose(nn2, i),
          exp(log(1-p) * (nn2-i) + log(p) * i) * (choose(nn2, i) - intf.n[i+1]),
          exp(log(1-p) * (ws-i) + log(p) * i) * choose(ws, i))
    })
}

## Iniitally defined only for heterodimers; Assumption is
## that there are size valid positions that can be a dincleotide
## the expected frequency of the dinucleotide is prob
## nucs and nucs.f may be used for homodimers where the
## independence of events is not clear (if AA, then probability
## of dn at -1 and +1 is simply the frequency of A; 
## this may disturb the distribution as probabilities increase
## with increased occurences of the dinucleotide.
dimer.dbinom <- function(size, prob, nucs=NULL, nucs.f=NULL){
    ## we can select up to size/2 heterodinucleotides, but we will consider
    ## the full set of size for the sake of completion.
    ##
    ## When one site is selected it decreases the number of available sites
    ## by 2 or 3, depending on whether it is odd, even, or terminal (the last
    ## available site). Hence, the number of valide sites is:
    ## N - (2*n.o + 3*n.e - e) 
    ##
    ## where e is is 2 if size is even (i.e. the terminal one reduces by -1 rather than -3
    ## if odd it reduces by -1 rather than -2
    N <- size
    e <- ifelse(N %% 2 == 0, 2, 1)
    ## n.o and n.e, can generally be assumed to be evenly selected
    ## if n is even; otherwise, odd is one larger.
    ## n.o <- ceiling( n / 2 )
    ## n.e <- floor( n/ 2)
    ## m is the number of available sites from which k will be selected
    m <- N - (2*n.o + 3*n.e - e/N)
    ## but this reduces too quickly. In fact we need to consider the full distribution
    ## of number of even and odd sites. sr, the reduction in sites:
    N.o <- ceiling(N / 2)
    N.e <- floor(N / 2)
    sr <- sapply(1:N, function(n){
        n.o <- ceiling(n/2)
        n.e <- n-n.o
        ## though I would prefer the equation to throw this up naturally
        if(n > N.o)
            return(rep(0, 6))
        p.odd.n <- dhyper(0:n, N.o, N.e, n)
        ## This gives me the number probability of selecting 0:n odd sites
        ## I will also want the probablity of k odd sites interfering with n-k even sites
        ## The sum of such probabilities will give me the probability of selections being
        ## illegal.
        p.allowed <- dhyper(0, 0:n , N.e-0:n, n:0) * p.odd.n
        ## Each odd position has a single even position that it interferes with
        ## Hence the probability of a given number of odd positions 
        ## The total number of ways of selecting n from N
        n.sel <- choose(N, n)
##        p.binomial <- n.sel * exp( log(prob) * n + log(1-prob) * (N-n) )
        c(n, n.o, n.e,
          n.sel * exp( log(prob) * n + log(1-prob) * (N-n) + log(sum(p.allowed))),
          binc=n.sel, p.al=sum(p.allowed))
    })
    ## add the first column for none selected.. 
    tmp=cbind(c(n=0, n.o=0, n.e=0, exp(log(1-prob) * n), 1, 1.0 ),
              sr)
    rownames(tmp) <- c('n', 'n.o', 'n.e', 'p', 'binc', 'p.al')
    tmp
}

## more simply..
matrix.mod <- function(ws, dn.p, homodimer=FALSE, fudge=1){
    ws <- ws + 1
##    dn.p <- dn.exp['TA']
    odd.n <- ceiling(ws/2)
    even.n <- floor(ws/2)
    bin.odd <- dbinom( 0:odd.n, odd.n, prob=dn.p )
    bin.even <- t(sapply( 0:odd.n, function(n.o){
        n.o <- ifelse(n.o == 0, 0, n.o-1)
        if(n.o >= even.n)
            return(c(1, rep(0, even.n)))
        f <- n.o / odd.n ## the frequency of sites blocked.
        left.n <- n.o * f ## the expected number of sites with a left neighbour
        block.n <- round(n.o - left.n)
        ## proportion that have no neighbour to their left
        ## we include the 0th position as a boundary condition
##        decrement <- round(n.o + n.o * f)
        decrement <- round((block.n + n.o) * fudge)
        decrement <- n.o
        n <- ifelse(homodimer, n, even.n - decrement)
        if(n < 1)
            return(c(1, rep(0, even.n)))
        c(dbinom( 0:n, n, prob=dn.p ), rep(0, even.n-n))
    }))
    p <- bin.odd * bin.even ## this has a length of 101, not 200.
    p.on <- row(p) - 1
    p.ev <- col(p) - 1
    p <- tapply(p, p.on + p.ev, sum)
    p
}
## this now produces perfect binomial distributions.. !! what?? 

## Another attempt. 
## I noted by considering the graph of all possible
## sequence paths (that has 4^n ends) that:
##
## a[i+1] <- (a[i] + b[i])/4
## b[i+1] <- (b[i] * 4) - (a[i] + b[i])/4
## i.e.
## b[i+1] <- (b[i] * 4) - a[i+1]
##
## Where a and b are the number of paths that terminate
## with a match (a), or have no matches at all for
## a given length of nucleotides.
## And the dinucleotide is a heterodimer. 
## a[1] <- 1
## b[1] <- 15
## if nucleotide composition is homogeneous. If not
## the equations become:
## a[1] <- f[1] * f[2] * 16
## b[1] <- 16 - a[1]
## and
## a[i+1] <- (a[i] + b[i]) * 4 * f[1] * f[2]
## b[i+1] <- b[i] * 4 - a[i+1]
## 
## If the dimer is a homodimer, then, the following is true
## a[1] <- 1
## b[1] <- 15
## a[i+1] <- (a[i] + b[i]) * 3 * 4 / 16
## b[i+1] <- (b[i] * 4) - a[i+1]
## This is a slight reduction; in fact it should not decrease more
## but be constant. It would be better to write a single expression
## that can be used for all:
## the probability that d[j] is a match is:
## f1[i] * f2[j]
## where f1 and f2 are the frequencies of nucleotide 1 and 2
## at positions i and j.
## if the last position was a match, then
## f1[i] is (k.nucs[1] == k.nucs[2]) i.e. 0 or 1
kp1.mod <- function(mnf, w, k.nucs){
    if(length(k.nucs) != 2)
        stop("only dimers working at the moment")
    n.f <- mnf[ k.nucs ]
    ## two matrices for position i and j (j=1+1)
    ## each matrix has two columns; the first one gives the
    ## counts of paths that end in a non-match (old)
    ## the second gives the ones that end in a match (new).
    ## the counts for 0 matches are given by the equation for
    ## b[i] above. Here the transition is from b[i] -> b[j]
    m.i <- matrix(0, nrow=w, ncol=2)
    colnames(m.i) <- c("old", "new")
    m.j <- m.i
    m.i[1,"old"] <- 16 * (1 - n.f[1] * n.f[2])
    m.i[2,"new"] <- 16 * n.f[1] * n.f[2]
    a <- rep(0, w)
    b <- a
    a[1] <- m.i[2,2]
    b[1] <- m.i[1,1]
    ## then we define m.j for the rest of the window;
    for(j in 2:w){
        i <- j-1
        a[j] <- (a[i] + b[i]) * 4 * n.f[1] * n.f[2]
        b[j] <- b[i] * 4 - a[j]
        m.j[1,"old"] <- b[j]
        m.j[2,"new"] <- a[j]
        ## then we transfer the counts for all chains that end in
        ## a match. We multiply by 4; this may not be necessary, but
        ## I don't know if the recursive rule works but it looks like
        ## I can remove the multiplication since it is done for both a and b
        m.j[,"old"] <- m.j[,"old"] + m.i[,"new"] * 4
        ## chains with matches, but which end in mismatches can then
        ## be transferred. We can probably consider them all equivalent and simply do:
        ## f[1] * f[2] * the numbers.
        ## note that I already set the m.j[2,"new"] above. 
        ## from the value of a[j]
        ## If I'm strict about it, I should use the recursive function for all
        ## "old" -> "new" transitions. But this would require holding the age of
        ## chains terminating in a match. This is not super difficult, but would
        ## be easier to implement in C with proper data structures.
        ## chain.
        ## But in general we should not simply use n.f[1] * n.f[2] because
        ## we know that the previous dinucleotide was not k.nucs[1:2]
        ## And this means that the probability is higher; at least a[2] / (a[2] + b[2])
        ## Ideally we would record the number of non-matches for every chain.
        new.n <- m.i[3:w-1,"old"] * ( a[2] / (a[2] + b[2])) ## n.f[1] * n.f[2]
        old.n <- m.i[3:w-1,"old"] - new.n
        m.j[3:w,"new"] <- m.j[3:w,"new"] + new.n * 4
        m.j[3:w-1,"old"] <- m.j[3:w-1,"old"] + old.n * 4
        ## and I think I can switch m.i and m;
        m.i <- m.j
        m.j[,] <- 0
    }
    list(m=m.i, a=a, b=b)
}

## For a dinucleotide (dn) composed of nucleotides n1 and n2
## that are present at frequencies f1 and f2.
## For a naive sampling the frequency of dn
## is f1 * f2
## The non-dn fraction is 1 - (f1 * f2)
## We wish to define the frequency of n1 in the non-dn fraction
## This is different if n1 == n2 or not.
## if n1 != n2 then it is simply
## f1 / (1 - f1 * f2)
## since none of (f1 * f2) can end in n1, but the total
## fraction of n1 in position 2 must be f1.
##
## if n1 == n2
## (f1 - f1 * f2) / (1 - f1 * f2)
##
## We can refer to this as f1.h; (f-hat)
## events since the last positive dn.
## a[i] is the proportion of positive events.
## a[i] is simply f.h * f2, as f2 is independently sampled.
## b[i] is simply 1 - a[i]
##
## And f1.h is defined from the a[i-1]
## 
## We then define a (the proportion of matching dinucleotides)
## and b for each step.
## mnf is a named vector of nucleotide frequencies
## w is the window width
## k.nucs are the nucleotides of the dimer
## freq : if true give frequencies; if false give counts
ab.chain <- function(mnf, w, k.nucs){
    a <- rep(0, w)
    b <- a
    ## m is a multipler; if freq is false it is 4, otherwise
    ## 1. It defines whether we report counts or frequencies
    ## h, defines if this is a homodimer or not
    h  <- as.numeric( k.nucs[1] == k.nucs[2] )
    f1 <- mnf[k.nucs[1]]
    f2 <- mnf[k.nucs[2]]
    a[1] <- f1 * f2 * 16
    b[1] <- 16 - a[1] ## (1 - f1 * f2)
    ## then we define m.j for the rest of the window;
    for(j in 2:w){
        i <- j-1
        ## 
        a[j] <- 4 * (a[i] + b[i]) * (f1 - (h * a[i]/(a[i]+b[i]))) * f2
        b[j]  <- 4 * b[i] - a[j]
    }
    list(a=a, b=b)
}


ls
## This is a simplified version; it uses the same logic as above,
## but does away with the complicated a and b chains. These were
## useful for me to define the problem, but I think they may not be
## needed.
## the probability that d[j] is a match is:
## f1[i] * f2[j]
## where f1 and f2 are the frequencies of nucleotide 1 and 2
## at positions i and j.
## note that f2[j] is always just f2 as it is randomly chosen.
## f1[i] depends on whether the last dinucleotide matched.
## 
## if the last position was a match, then
## f1[i] is (k.nucs[1] == k.nucs[2]) i.e. 0 or 1
## so the probability of a match is just:
## (k.nucs[1] == k.nucs[2]) * f2
##
## If however the position was not a match, then f1[i] is
## f1[i] <- f1 / (1 - f1 * f2)
## or if homodimer the probability is decreased
## f1[i] <- (1-f1) * f2 / (1-f1 * f2)
## we can use these equations in a similar manner.
kp2.mod <- function(mnf, w, k.nucs, freq=TRUE){
    if(length(k.nucs) != 2)
        stop("only dimers working at the moment")
    n.f <- mnf[ k.nucs ]
    t.m <- ifelse(freq, 1, 4)
    ## two matrices for position i and j (j=1+1)
    ## each matrix has two columns; the first one gives the
    ## counts of paths that end in a non-match (old)
    ## the second gives the ones that end in a match (new).
    ## the counts for 0 matches are given by the equation for
    nr <- w+1
    m.i <- matrix(0, nrow=nr, ncol=2)
    colnames(m.i) <- c("old", "new")
    m.i[1,"old"] <- t.m^2 * (1 - n.f[1] * n.f[2])
    m.i[2,"new"] <- t.m^2 * n.f[1] * n.f[2]
    homo = as.numeric(k.nucs[1] == k.nucs[2])
    ## then we define m.j for the rest of the window;
    for(j in 2:w){
        ## new column: these end in a match. 
        ## if homodimers then promotion can happen
        new.n1 <- m.i[,'new'] * homo * n.f[2]
        old.n1 <- m.i[,'new'] - new.n1
        ## old column; these do not end in match
        new.n2 <- m.i[,'old'] * n.f[2] * (n.f[1] - homo * n.f[1] * n.f[2]) / (1-n.f[1] * n.f[2])
        old.n2 <- m.i[,'old'] - new.n2
        m.i[2:nr,'new'] <- t.m * (new.n1[-nr] + new.n2[-nr])
        m.i[,'old'] <- t.m * (old.n1 + old.n2)
    }
    m.i
}
### kp2.mod almost recreates the random results; however it
### does not completely model the homodimers. I think that I
### do actually need to use the ab chain concept and keep a matrix
### with several columns denoting the ages of non-matches
### and then make use of the ab-chain to modify the resulting
### statistic.
### It seems above that the penalty for homodimers is too high after
#### the first step. And this has an effect on the numbers of 
### low counts.

### I found that the difference between the inferreed and observed
### distributions were worse for short windows.
### That seems to have been due to a change to the code that made things worse.
### Hence lets try different windows, but with a shorter random sequence so that it
### is a bit quicker.
### also include a homogeneous all 1/4 set of values.
rnd.1 <- substring(lp.rnd, 1, 5e6)
rnd.2 <- paste0( sample(kmer.seq(1), 5e6, replace=TRUE, prob=rep(1/4, 4)), collapse="")
rnd.2.2 <- paste0( sample(kmer.seq(1), 5e6, replace=TRUE, prob=rep(1/4, 4)), collapse="")
rnd.2.3 <- paste0( sample(kmer.seq(1), 5e6, replace=TRUE, prob=rep(1/4, 4)), collapse="")

## count two k-mers
rnd.1.wc.2 <- window.kmer.dist( rnd.1, dinucleotides, 2, freq=TRUE )
rnd.2.wc.2 <- window.kmer.dist( rnd.2, dinucleotides, 2, freq=TRUE )
## do for the very long sequence. (~1 Gbp)
rnd.l.wc.2 <- window.kmer.dist( lp.rnd, dinucleotides, 2, freq=TRUE )

rnd.1.exp <- lapply( strsplit(dinucleotides, ""), function(x){
    kp2.mod( lp.mnf$f, 2, x, freq=TRUE)})
names(rnd.1.exp) <- dinucleotides

even.f <- lp.mnf$f
even.f[] <- 1/4
rnd.2.exp <- lapply( strsplit(dinucleotides, ""), function(x){
    kp2.mod( even.f, 2, x, freq=TRUE)})
names(rnd.2.exp) <- dinucleotides

## to compare we can do:
homodi <- sapply(strsplit(dinucleotides, ""), function(x){ x[1] == x[2] })
pch <- ifelse(homodi, 19, 1)

par(mfrow=c(3,1))
plot(col(rnd.1.wc.2$dist), log10(sapply(rnd.1.exp, rowSums) / rnd.1.wc.2$dist ), col=row(rnd.1.wc.2$dist), pch=pch[col(rnd.1.wc.2$dist)]
   , ylab='log ratio', xlab="kmer")
abline(h=0, lty=2)
##
##
plot(col(rnd.l.wc.2$dist), log10(sapply(rnd.1.exp, rowSums) / rnd.l.wc.2$dist ), col=row(rnd.l.wc.2$dist), pch=pch[col(rnd.l.wc.2$dist)]
   , ylab='log ratio', xlab="kmer")
abline(h=0, lty=2)
##
plot(col(rnd.l.wc.2$dist), log10(rnd.1.wc.2$dist / rnd.l.wc.2$dist ), col=row(rnd.l.wc.2$dist), pch=pch[col(rnd.l.wc.2$dist)]
   , ylab='log ratio', xlab="kmer")
abline(h=0, lty=2)

## that actually says that we are pretty much correct here. The longer (much) random sequence
## gives better estimates. Now lets try for a window with 3 k-mers and see if we get what we expect.
rnd.1.wc.3 <- window.kmer.dist( rnd.1, dinucleotides, 3, freq=TRUE )
rnd.2.wc.3 <- window.kmer.dist( rnd.2, dinucleotides, 3, freq=TRUE )
##
rnd.2.2.wc.3 <- window.kmer.dist( rnd.2.2, dinucleotides, 3, freq=TRUE )
rnd.2.3.wc.3 <- window.kmer.dist( rnd.2.3, dinucleotides, 3, freq=TRUE )

## do for the very long sequence. (~1 Gbp)
rnd.l.wc.3 <- window.kmer.dist( lp.rnd, dinucleotides, 3, freq=TRUE )

rnd.1.exp.3 <- lapply( strsplit(dinucleotides, ""), function(x){
    kp2.mod( lp.mnf$f, 3, x, freq=TRUE)})
names(rnd.1.exp.3) <- dinucleotides

rnd.2.exp.3 <- lapply( strsplit(dinucleotides, ""), function(x){
    kp2.mod( even.f, 3, x, freq=TRUE)})
names(rnd.2.exp.3) <- dinucleotides

par(mfrow=c(3,1))
plot(col(rnd.1.wc.3$dist), log10(sapply(rnd.1.exp.3, rowSums) / rnd.1.wc.3$dist ), col=row(rnd.1.wc.3$dist), pch=pch[col(rnd.1.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)
##
plot(col(rnd.l.wc.3$dist), log10(sapply(rnd.1.exp.3, rowSums) / rnd.l.wc.3$dist ), col=row(rnd.1.wc.3$dist), pch=pch[col(rnd.1.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)
plot(col(rnd.l.wc.3$dist), log10(rnd.1.wc.3$dist / rnd.l.wc.3$dist ), col=row(rnd.1.wc.3$dist), pch=pch[col(rnd.1.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)

par(mfrow=c(3,1))
plot(col(rnd.2.wc.3$dist), log10(sapply(rnd.2.exp.3, rowSums) / rnd.2.wc.3$dist ), col=row(rnd.1.wc.3$dist), pch=pch[col(rnd.1.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)
##
plot(col(rnd.2.2.wc.3$dist), log10(sapply(rnd.2.exp.3, rowSums) / rnd.2.2.wc.3$dist ), col=row(rnd.2.wc.3$dist), pch=pch[col(rnd.2.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)
plot(col(rnd.2.3.wc.3$dist), log10(sapply(rnd.2.exp.3, rowSums) / rnd.2.3.wc.3$dist ), col=row(rnd.2.wc.3$dist), pch=pch[col(rnd.2.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)
###
### here we see a systematic difference:
### for homodimers:
#### expected counts of 1 (red) are lower than they should be. Counts of 2 are higher, counts of 3 correct
### for heterodimers:
### the same patter is observed, but the difference for count of 1 (red) is much less (and no count 3)

### First confirm. Are we getting the correct counts?
s.t <- c("AAAATGAAATTAA",    ## 3: 1, 2: 3, 1: 4 0: 2
         "AAAATGAAATTAAA",   ## 3: 1, 2: 4, 1: 4 0: 2
         "AAAATGAAATTAAT",   ## 3: 1, 2: 3, 1: 5 0: 2
         "TAAAATGAAATTAAT")   ## 3: 1, 2: 4, 1: 5 0: 2
tmp <- lapply(s.t, window.kmer.dist, kmers="AA", window=3, freq=FALSE, ret.flag=1L)
## this does seem to give the correct numbers for all:
## Then we try to work out what we should get...
## start with simple situation of even.f

rnd.2.wc.5 <- window.kmer.dist( rnd.2, dinucleotides, 5, freq=TRUE )
##
rnd.2.2.wc.5 <- window.kmer.dist( rnd.2.2, dinucleotides, 5, freq=TRUE )
rnd.2.3.wc.5 <- window.kmer.dist( rnd.2.3, dinucleotides, 5, freq=TRUE )

rnd.2.wc.4 <- window.kmer.dist( rnd.2, dinucleotides, 4, freq=TRUE )
##
rnd.2.2.wc.4 <- window.kmer.dist( rnd.2.2, dinucleotides, 4, freq=TRUE )
rnd.2.3.wc.4 <- window.kmer.dist( rnd.2.3, dinucleotides, 4, freq=TRUE )

## window of 7
rnd.2.wc.7 <- window.kmer.dist( rnd.2, dinucleotides, 7 )

par(mfrow=c(3,1))
plot(col(rnd.2.wc.3$dist), log2(sapply(rnd.2.exp.3, rowSums) / rnd.2.wc.3$dist ), col=row(rnd.2.wc.3$dist), pch=pch[col(rnd.2.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)

##
##
plot(col(rnd.l.wc.3$dist), log10(sapply(rnd.2.exp.3, rowSums) / rnd.l.wc.3$dist ), col=row(rnd.l.wc.3$dist), pch=pch[col(rnd.l.wc.3$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)
##
plot(col(rnd.l.wc.3$dist), log10(rnd.2.wc.3$dist / rnd.l.wc.3$dist ), col=row(rnd.l.wc.3$dist), pch=pch[col(rnd.l.wc.2$dist)]
   , ylab='log ratio', xlab="kmer", cex=1.5)
abline(h=0, lty=2)



i <- i:4
par(mfrow=c(1,1))
plot(i-1, log10(rnd.1.wc.3$dist[,'AA']), cex=1.5, ylim=range(log10(rnd.1.wc.3$dist), finite=TRUE))
points(i-1, log10(rnd.l.wc.3$dist[,'AA']), cex=1, pch=19)
lines(i-1, log10(rowSums(rnd.1.exp.3[['AA']])))
##
points(i-1, log10(rnd.1.wc.3$dist[,'CC']), cex=1.5, col=2)
points(i-1, log10(rnd.l.wc.3$dist[,'CC']), cex=1, pch=19, col=2)
lines(i-1, log10(rowSums(rnd.1.exp.3[['CC']])), col=2)
##
points(i-1, log10(rnd.1.wc.3$dist[,'TA']), cex=1.5, col=3)
points(i-1, log10(rnd.l.wc.3$dist[,'TA']), cex=1, pch=19, col=3)
lines(i-1, log10(rowSums(rnd.1.exp.3[['TA']])), col=3)
##
points(i-1, log10(rnd.1.wc.3$dist[,'CG']), cex=1.5, col=4)
points(i-1, log10(rnd.l.wc.3$dist[,'CG']), cex=1, pch=19, col=4)
lines(i-1, log10(rowSums(rnd.1.exp.3[['CG']])), col=4)

#### After extensive testing and manual calculations it is apparent
#### that we do need to consider the change in transition probabilities
#### To investigate this, I will use a recursive function that extracts the
#### selection graph with associated transfer probabilities;

## This is only for dinucleotides
transition.p <- function(mnf, nucs, depth, freq=FALSE){
    if(length(nucs) != 2)
        stop("Function only defined for dinucleotides")
    ## n1.f, should be the frequency of the first nucleotide of the
    ## the dimer in the terminal position; the frequency of dimers
    ## is then n1.f * mnf[ nucs[2] ]
    ## we will also define f1 and f2 for convenience
    ## and mlt which mulitplies if we want to look at total numbers
    ## of nodes.
    ## x and y give the depth and a y position of the node.
    descend.tree <- function(n, n1.f, x, y, match.n, age){
        ## f1 and f2 are defined in the calling scope
        ## as are n1 and n2 the two nucleotides of the dimer
        ## mlt is 1 if freq == FALSE, otherwise 4
##        browser()
        if(x == 0){ ## the first step: initialise using just f1 and f2
            dn.true <- f1 * f2 * n * mlt
            dn.false <- n * mlt - dn.true
            n1.f <- f1
        }else{
            dn.true <- n1.f * f2 * mlt * n
            dn.false <- n * mlt - dn.true
        }
        ## modify an object in the calling scope called graph:
        graph <<- rbind(graph, c(n * mlt, dn.true, dn.false, x, y, 0, 0, match.n, age))
        cr <- nrow(graph)
        if(x == depth)
            return(cr)
        ## then work out what is the frequency of n1 in the two
        ## branches. This depends on whether the nucleotide is
        ## a homodimer or heterodimer
        n <- n * mlt
        if(n == 0){
            return(cr)
        }
        true.n1.f <- ifelse(n1 == n2, 1, 0)
        false.n1.f <- ifelse(n1 == n2, (n * f1 - dn.true)/dn.false, f1 / (1 - dn.true/n))
        graph[cr, 'down'] <<- descend.tree( dn.false, false.n1.f, x + 1, y - 1/2^x, match.n, age+1 )
        graph[cr, 'up'] <<- descend.tree( dn.true, true.n1.f, x + 1, y + 1/2^x, match.n+1, 0 )
        return(cr)
    }
    n1 <- nucs[1]; n2 <- nucs[2];
    f1 <- mnf[n1]; f2 <- mnf[n2];
    mlt <- ifelse(freq, 1, 4)
    graph <- matrix(ncol=9, nrow=0)
    colnames(graph) <- c("n", "m", "mm", "x", "y", "up", "down", "match.n", "age")
    descend.tree( mlt, 0, 0, 0, 0, 2 )
    as.data.frame(graph)
}

draw.graph <- function(graph, cex=1, adj.y=0.3){
    with(graph, plot(x, y, type='n'))
    with(graph, text(x, y, round(m), cex=cex, adj=c(0.5, 0-adj.y)))
    with(graph, text(x, y, round(mm), cex=cex, adj=c(0.5, 1 + adj.y)))
    sep.l <- with(graph, pmax( strwidth(round(m), cex=cex), strwidth(round(mm), cex=cex))) / 2
    with(graph, segments(x - sep.l, y, x + sep.l, y ))
    l.h <- strheight("A", cex=cex) / 1.5
    lw.i <- strwidth("M", units="inches") ### to adjust arrows
    b <- graph$up > 0
    ar.sep <- max(sep.l) * 1.2
    with(graph, arrows(x[b]+ar.sep, y[b]+l.h,
                     x[up[b]]-ar.sep, y[up[b]],
                     length=lw.i, angle=20))
    with(graph, arrows(x[b]+ar.sep, y[b]-l.h,
                     x[down[b]] - ar.sep, y[down[b]],
                     length=lw.i, angle=20))

}
    

tmp1 <- transition.p( even.f, c("T", "A"), 7, freq=FALSE )
tmp2 <- transition.p( even.f, c("A", "A"), 7, freq=FALSE )

pdf("transition_graphs.pdf", width=12, height=50)
par(mfrow=c(2,1))
draw.graph(tmp1, cex=0.5, adj.y=0.3)
draw.graph(tmp2, cex=0.5, adj.y=0.3)
dev.off()

### and check how it compares to the sequences that we normally get.
graph.to.dist <- function(graph, x){
    graph <- graph[ graph$x == x, ]
    counts <- with(graph, rbind( cbind(match.n, mm), cbind(match.n+1, m) ))
    tapply( c(0:x, counts[,2]), c(rep(0,x+1), counts[,1]), sum )
}

## turns out that the 7th position gives me one more, so a count of up to 8
rnd.2.wc.8 <- window.kmer.dist( rnd.2, dinucleotides, 8, freq=TRUE )
rnd.2.2.wc.8 <- window.kmer.dist( rnd.2.2, dinucleotides, 8, freq=TRUE )
rnd.2.3.wc.8 <- window.kmer.dist( rnd.2.3, dinucleotides, 8, freq=TRUE )

tmp1.c <- graph.to.dist(tmp1, 7)
tmp2.c <- graph.to.dist(tmp2, 7)

rnd.1.wc.8 <- window.kmer.dist( rnd.1, dinucleotides, 8, freq=TRUE )
rnd.l.wc.8 <- window.kmer.dist( lp.rnd, dinucleotides, 8, freq=TRUE )

tmp3 <- transition.p( lp.mnf$f, c("T", "A"), 7, freq=FALSE )
tmp4 <- transition.p( lp.mnf$f, c("A", "A"), 7, freq=FALSE )
tmp3.c <- graph.to.dist(tmp3, 7)
tmp4.c <- graph.to.dist(tmp4, 7)

## Up to count of 8 seems to be pretty consistent. So far, so good.
## try for a longer window.
system.time(
    tmp5 <- transition.p( lp.mnf$f, c("T", "A"), 20, freq=FALSE )
)
##   user  system elapsed 
## 17.367   1.264  18.631 

## This is way, way, slower; as fewer branches end in 0xo
tmp6 <- transition.p( lp.mnf$f, c("A", "A"), 20, freq=FALSE )

tmp5.c <- graph.to.dist(tmp5, 20)
tmp6.c <- graph.to.dist(tmp6, 20)

rnd.1.wc.21 <- window.kmer.dist( lp.rnd, dinucleotides, 21, freq=TRUE )
## It still looks like this might be giving me the correct data sets.

## Have a look for some patterns that we might be able to use:
with(tmp5, tapply(1:length(n), x, function(i){ table( m[i] / n[i] ) }))
with(tmp6, tapply(1:length(n), x, function(i){ table( m[i] / n[i] ) }))

with(tmp5, tapply(1:length(n), x, function(i){ table( m[i] / n[i], age[i] ) }))
with(tmp6, tapply(1:length(n), x, function(i){ table( m[i] / n[i], age[i] ) }))

## this solves the complete problem:
kp3.mod <- function(mnf, w, k.nucs){
    if(length(k.nucs) != 2)
        stop("only dimers working at the moment")
    n1 <- k.nucs[1]; n2 <- k.nucs[2]
    f1 <- mnf[ n1 ]
    f2 <- mnf[ n2 ]
    ## we consider a set of probabilites of the next nucleotide resulting
    ## in a window. 
    ## Where we have:
    ## mp.init : The initial matching probablity. This starts at
    ##           f1 * f2, but then changes with each step. This is the
    ##           probability of leaving the 0-match branch. Only a single
    ##           value is required, but it will be updated each step
    ## mm.p[i] : A vector of values giving the probality of a mismatch of
    ##           age i counted from the last match. There can w-1 different
    ##           values. The relationship between the current mp.init and
    ##           these values depend on the whether the dinucleotide repeats or not
    ## m2m.p   : the probability of a consecutive match. This is either 0 or f2
    ##           depending on whether n1 == n2.
    mp.init <- f1 * f2
    tr.p <- rep(0, length=w+1)
    m2m.p <- ifelse(n1 == n2, f2, 0)
    ## we then use the update rules as set in the transition.p function defined above to
    ## explore the trees made in this case.
    m0.f <- 1  ## this is the initial frequency.
    ## we have a table of frequencies for different counts at different ages
    ## the rows indicate ages; we have an extra row and column
    counts <- matrix(0, nrow=w+1, ncol=w+1)
    counts.tmp <- counts
    ## rows indicate ages..
    tr.p[1] <- m2m.p ## this should be remove later: this is only for testing purposes
    for(i in 1:w){
##        tr.p[i+1] <- mp.init
        m.new <- mp.init * m0.f ## the fraction of new matches
        mm.f <- m0.f - m.new
        m0.f <- m0.f - m.new
        ## n1.f is the frequency of n1 in the second position in case of a mismatch
        ## for heterodimers this is higher than f2, for homodimers it is lower
        tr.pp <- tr.p[i]
        n1.f.1 <- ifelse( n1 == n2, (f1 - mp.init)/(1-mp.init), f1 / (1-mp.init) )
        n1.f.2 <- ifelse( n1 == n2, (f1 - tr.pp)/(1-tr.pp), f1 / (1-tr.pp) )
        mp.init <- n1.f.1 * f2
        tr.p[i+1] <- n1.f.2 * f2
        ## a temporary counts matrix:
        counts.tmp[1,1:i+1] <- colSums(counts[1:i,1:i,drop=FALSE] * tr.p[1:i])
        counts.tmp[1:i+1,1:i] <- counts[1:i,1:i] * (1-tr.p[1:i])
        counts.tmp[1,1] <- counts.tmp[1,1] + m.new
        counts <- counts.tmp
        counts.tmp[,] <- 0
    }
    list( dist=c(m0.f, colSums(counts)[-(w+1)]), tr=tr.p )
}

tmp5.t <- kp3.mod( lp.mnf$f, w=21, k.nucs=c("T", "A") )
tmp6.t <- kp3.mod( lp.mnf$f, w=21, k.nucs=c("A", "A") )

plot(tmp5.t$dist[1:13], tmp5.c/sum(tmp5.c))
abline(0,1) ## looks good

## but we have a potential problem here:
plot(tmp6.t$dist, tmp6.c/sum(tmp6.c))
abline(0, 1) ## looks good

plot(tmp6.t$dist, type='l')
points(1:length(tmp6.c), tmp6.c / sum(tmp6.c), col='red', type='l')
points(1:length( rnd.1.wc.21$dist[,'AA']), rnd.1.wc.21$dist[,'AA'], type='l', col='blue' )
## That looks perfect!

## long long .. making the sequence this way is extremely slow...
lp.rnd.ll <- paste(sample(kmer.seq(1), 2e9, replace=TRUE, prob=lp.mnf$f), collapse="")
lp.rnd.ll.wc <- window.kmer.dist( lp.rnd.ll, dinucleotides, 200 )
                    
## We can now try to do this for a longer window and see if I have cracked the nut..
tmp.7 <- kp3.mod( lp.mnf$f, w=200, k.nucs=c("T", "A"))
tmp.8 <- kp3.mod( lp.mnf$f, w=200, k.nucs=c("A", "A"))

pdf("dimer_random_distributions_TA_AA.pdf", width=7, height=7)
i <- 1:40
plot(i-1, tmp.7$dist[i], type='l')
points(i-1, lp.rnd.wc$dist[i,'TA']) ## ok, looking good
points(i-1, lp.rnd.ll.wc$dist[i,'TA'], pch=19, cex=0.8) ## ok, looking good
lines(i-1, tmp.8$dist[i], type='l', col='red')
points(i-1, lp.rnd.wc$dist[i,'AA'], col='red')
points(i-1, lp.rnd.ll.wc$dist[i,'AA'], col='red', pch=19, cex=0.8)
legend("topright", legend=c("TA", "AA", "calculated", "sim short", "sim long"),
       pch=c(NA, NA, NA, 1, 19), lwd=c(1, 1, 1, NA, NA), col=hsv(0, c(1,1,0.5,0.5,0.5), c(0, 1, 0, 0, 0)))
dev.off()

## lets do this for the full set of dinucleotides and then compare to the real distribution
## and the observed one.
expected.dimer.dists <- lapply( strsplit(dinucleotides, ""), function(x){
    kp3.mod( lp.mnf$f, w=200, k.nucs=x )
})
names(expected.dimer.dists) <- dinucleotides
## that's more or less instantaneous:

par(mfrow=c(4,4))
i <- 1:50
x <- i-1
for(dn in dinucleotides){
    dists <- cbind(lp.wc$dist[,dn], lp.rnd.ll.wc$dist[,dn], expected.dimer.dists[[dn]]$dist)
    plot(x, dists[i,1], type='l', ylim=range(dists[i,]), xlab="number of matches", ylab="frequency", main=dn, lwd=2)
    lines(x, dists[i,2], col='red', lwd=2)
    lines(x, dists[i,3], col=rgb(0, 0, 1, 1), lty=2)
}

ls
## lets check distributions of random numbers:
## this is not so fast, but should be sufficient to give us an idea
## actually for this, i should do 100, but lets not care
tmp.200 <- lapply(1:200, function(j){
    t(sapply(1:1000, function(i){
        diff(c(0, sort(sample(1:200, j)), 201))
    }))
})

tmp.100 <- lapply(1:100, function(j){
    t(sapply(1:1000, function(i){
        diff(c(0, sort(sample(1:100, j)), 101))
    }))
})

## get the distribution of l == 1 and right d == 1 for all rows
tmp.200.1 <- lapply(tmp.200, function(x){
    j <- 2:ncol(x)
    i <- j-1
    b <- x == 1
    cbind( left=rowSums(b[,i,drop=FALSE]), right=rowSums(b[,j,drop=FALSE]), both=rowSums(b[,i,drop=FALSE] & b[,j,drop=FALSE]),
          single=rowSums(!b[,i,drop=FALSE] & !b[,j,drop=FALSE]),
          left.only=rowSums(b[,i,drop=FALSE] & !b[,j,drop=FALSE]), right.only=rowSums(!b[,i,drop=FALSE] & b[,j,drop=FALSE]))
})

tmp.100.1 <- lapply(tmp.100, function(x){
    j <- 2:ncol(x)
    i <- j-1
    b <- x == 1
    cbind( left=rowSums(b[,i,drop=FALSE]), right=rowSums(b[,j,drop=FALSE]), both=rowSums(b[,i,drop=FALSE] & b[,j,drop=FALSE]),
          single=rowSums(!b[,i,drop=FALSE] & !b[,j,drop=FALSE]),
          left.only=rowSums(b[,i,drop=FALSE] & !b[,j,drop=FALSE]), right.only=rowSums(!b[,i,drop=FALSE] & b[,j,drop=FALSE]))
})

tmp.200.1.f <- as.data.frame(t(sapply(tmp.200.1, function(x){
    colSums(x) / sum(x[,3:6])
})))

tmp.100.1.f <- as.data.frame(t(sapply(tmp.100.1, function(x){
    colSums(x) / sum(x[,3:6])
})))

f <- (1:100 + 2)/102
with(tmp.100.1.f, plot(1:length(left), single))
f <- (1:100 + 2)/102
with(tmp.100.1.f, lines(1:length(left), left))
with(tmp.100.1.f, lines(1:length(left), right))
points(1:length(f), (1-f^2), cex=0.5, col='red')

x <- 1:100
tmp.100.intcpt <- t(apply( tmp.100.1.f, 2, function(y){
    summary(lm( y ~ x ))$coefficients[1,] ## 5.000e-03
}))
tmp.100.x <- t(apply( tmp.100.1.f, 2, function(y){
    summary(lm( y ~ x ))$coefficients[2,] ## 5.000e-03
}))

with(tmp.100.1.f, plot(1:length(left), left, type='l'))
with(tmp.100.1.f, lines(1:length(left), right, type='l'))
with(tmp.100.1.f, lines(1:length(left), both, type='l', col='red'))
with(tmp.100.1.f, lines(1:length(left), single, type='l', col='purple'))
with(tmp.100.1.f, lines(1:length(left), sqrt(single), type='l', col='purple'))

par(mfrow=c(3,1))
with(tmp.100.1.f, plot( x, log2( (x / 100) / left ) ) )
with(tmp.100.1.f, points( x, log2( (x / 100) / right ), col='red' ) )
abline(h=0)
with(tmp.100.1.f, plot( x, log2( (x / 101) / left ) ) )
with(tmp.100.1.f, points( x, log2( (x / 101) / right ), col='red' ) )
abline(h=0)
with(tmp.100.1.f, plot( x, log2( (x / 102) / left ) ) )
with(tmp.100.1.f, points( x, log2( (x / 102) / right ), col='red' ) )
abline(h=0)
### Interestingly, 1 / 100, gives the best result here. I guess that's the effective
### frequency, as the borders are fixed and can't move. It should suffice.

## counts for left, right and double? But left and right should be the same, but they are not?
## double is the probability of both left and right, it should be equal to left * right,
## single is neither left nor right. 1 - (1-left) * (1-right)
par(mfrow=c(1,1))
with(tmp.100.1.f, plot( x, single))
## 
lines( x, (1 - x/100) * (1 - x/100), col='red')
with(tmp.100.1.f, points( x, both ))
lines( x, (x/100)*(x/100), col='blue')
with(tmp.100.1.f, points( x, left))
lines(

ws <- 200
for(dn in names(dn.exp)){
    homodimer <- length(unlist(strsplit(dn, ""))) == 1
    dn.p <- matrix.mod( ws, dn.exp[dn], homodimer=homodimer, fudge=1)
    ## the normal pbinomial is not distinguishable..
    i <- 1:40
    x <- i-1
    dbm <- dbinom( 0:ws, ws, prob=dn.exp[dn] )
    par(mfrow=c(1,2))
    with(lp.rnd.wc, plot( x, dist[i,dn], type='l', main=dn,
                         ylim=range(c(dbm, dist[,dn], dn.p)), lwd=3))
    mx <- which.max(lp.rnd.wc$dist[,dn])
    abline(v=x[mx], lty=2)
    lines(x, dn.p[i], type='l', col='red', lwd=3)
    ## and the normal dbinomial
    lines(x, dbm[i], type='l', col='blue', lwd=3)
    ## 
    with(lp.rnd.wc, plot( x, cumsum(dist[i,dn]), type='l', main=dn,
                         ylim=range(c(0,1)), lwd=3))
    lines(x, cumsum(dn.p[i]), type='l', col='red', lwd=3)
    ## and the normal dbinomial
    lines(x, cumsum(dbm[i]), type='l', col='blue', lwd=3)
    inpt <- readline("next: ")
}


tmp <- t(dimer.dbinom( 200, dn.exp['AT'] ))

plot( lp.rnd.wc$dist[,'AT'] )
lines(1 + 1:nrow(tmp), (tmp[,4]/sum(tmp[,4])))
dn <- "AT"
n <- 200
bin.1 <- dbinom( 0:n, n, p=dn.exp[dn] )
lines(1:length(bin.1), bin.1, col='red')

binom.err <- lapply( colnames(lp.rnd.wc$dist), function(dn){
    bin <- dbinom( 0:n, n, p=dn.exp[dn] )
    p <- exp( log(dn.exp[dn]) * 0:n + log(1-dn.exp[dn]) * n:0 )
    bin.c <- choose(n, 0:n)
    data.frame(rnd=lp.rnd.wc$dist[,dn], bin=bin, r=lp.rnd.wc$dist[,dn]/bin,
          binc=bin.c, p=p, rpb=lp.rnd.wc$dist[,dn]/p, bpb=bin/p)
})
names(binom.err) <- colnames(lp.rnd.wc$dist)

for(dn in names(dn.exp)){
    tmp <- as.data.frame(t(dimer.dbinom( 200, dn.exp[dn] )))
##    par(mfrow=c(2,2))
    layout(cbind(1:2, 3:3))
    i <- 1:40
    with(binom.err[[dn]], plot(i-1, rnd[i], type='b', ylim=range(c(rnd, bin)), main=dn))
    with(binom.err[[dn]], lines(i-1, bin[i], type='l', col='red'))
    with(tmp[i,], lines(n, p, col='blue'))
    with(tmp[i,], lines(n, p/sum(p), col='blue', type='b', cex=0.5))
    with(binom.err[[dn]], plot(i-1, log2(1/r[i]), type='b', ylim=range(
                                                              log2(c(1/r[i], tmp$p[i] / rnd[i], tmp$p[i] / bin[i])), finite=TRUE)))
    abline(h=0, lty=3)
##    with(binom.err[[dn]], lines(i-1, log2(rpb[i]/(tmp$p.al[i]*tmp$binc[i])), type='l', col='red'))
    with(binom.err[[dn]], lines(i-1, log2(tmp$p[i]/bin[i]), type='l', col='red'))
    with(binom.err[[dn]], lines(i-1, log2(tmp$p[i]/rnd[i]), type='l', col='blue'))
    lines( i-1, log2(tmp$p.al[i]), col='gold', lwd=2, lty=2 )
    with(par(), plot.window(xlim=usr[1:2], ylim=c(0,1), xaxs='i'))
    lines( i-1, tmp$p.al[i], col='gold', lwd=2 )
    ## the following is promising; 
    with(binom.err[[dn]], plot(i-1, log2(rpb[i]), type='b', ylim=range( log2(c(rpb[i], tmp[i,'p.al'] * tmp[i,'binc'])), finite=TRUE)))
    with(tmp[i,], lines(n, log2(p.al * binc), col='red'))
    with(tmp[i,], lines(n, log2(binc), col='blue'))
##    with(binom.err[[dn]], plot(i-1, log2(rpb/bpb)[i], type='b'))
    ## and the following gives exactly the same plot as log2(r) as we
    ## expect..
    input <- readline("next: ")
}




## still does worse than normal binomila. And somehow seem to give values that
## are way to high. I'm missing some componnent.

ws <- nrow(lp.wc$dist) - 1
##par(mfrow=c(4,4))
par(mfrow=c(1,2))
par(pmar=c(5.1, 4.1, 4.1, 4.1))
for(dn in colnames(lp.wc$dist)){
    i <- 1:40
    nucs <- unlist(strsplit(dn, ""))
    homodimer <- length( table(nucs) ) == 1
    ##    n <- ifelse(homodimer, ws-1, ceiling((ws-1)/2))
    n <- ws
    bin.1 <- dbinom( 0:n, n, p=dn.exp[dn] )
    p <- model.window.2(ws, lp.mnf$f, nucs)
    ##    bin.2 <- dpoisbinom(0:ws, p)
    m5 <- mod5(dn.exp[dn], ws)
    bin.2 <- m5[5,]
    bin.3 <- dn.bin( n, dn.exp[dn], homodimer )
##    hyp <- dhyper( 0:n, n * dn.exp[dn], n * dn.exp[1-dn], n )
    ## bin.3 <- dbinom( 0:n, n, p=lp.rnd.dn.exp[dn] )
    ## pois <- dpois( 0:n, lambda=(n * lp.rnd.dn.exp[dn] ))
##    
    plot(i-1, lp.wc$dist[i,dn], main=dn, type='l', ylim=range(c(lp.wc$dist[,dn], lp.rnd.wc$dist[,dn], bin.1)), lwd=2 )
    lines(i-1, lp.rnd.wc$dist[i,dn], col='red', lwd=2)
##    lines(1:nrow(lp.rnd.2.wc$dist)-1, lp.rnd.2.wc$dist[,dn], col='red', lwd=2)
    lines(i-1, bin.1[i], col='blue', lwd=2)
    lines(i-1, bin.2[i], col='green', lwd=2)
    lines(bin.3[i,1], bin.3[i,'p'], col='purple', lwd=2)
    plot(i-1, cumsum(lp.wc$dist[i,dn]),
         main=dn, type='l', ylim=c(0,1), lwd=2 )
    lines(i-1, cumsum(lp.rnd.wc$dist[i,dn]), col='red', lwd=2)
##    lines(1:nrow(lp.rnd.2.wc$dist)-1, lp.rnd.2.wc$dist[,dn], col='red', lwd=2)
    lines(i-1, cumsum(bin.1[i]), col='blue', lwd=2)
    lines(i-1, cumsum(bin.2[i]), col='green', lwd=2)
    lines(bin.3[i,1], cumsum(bin.3[i,'p']), col='purple', lwd=2)
    delta.1 <- bin.1[i] - lp.rnd.wc$dist[i,dn]
    delta.2 <- cumsum(bin.1[i]) - cumsum(lp.rnd.wc$dist[i,dn])
    with(par(), plot.window(xlim=usr[1:2], ylim=range(c(delta.1, delta.2)), xaxs='i'))
    axis(4)
    lines(i, delta.1, col="blue", lty=2)
    lines(i, delta.2, col="red", lty=2)
    identify(1,1)
##    lines(bin.2[,1], bin.2[,2], col="green", lwd=2)
##    lines(1:nrow(lp.rnd.wc$dist)-1, bin.3, col="purple", lwd=2)
##    lines(0:n, pois, col="purple", lwd=2)
}


i <- 1:40
plot(i-1, lp.rnd.wc$dist[i,"AT"], col='black', lwd=2, type='b')
lines(i-1, lp.rnd.wc$dist[i,"TA"], col='red', lwd=2)
lines(i-1, lp.rnd.wc$dist[i,"TT"], col='green', lwd=2, type='b')
lines(i-1, lp.rnd.wc$dist[i,"AA"], col='blue', lwd=2)
m5 <- mod5(dn.exp["AT"], ws)
lines(i-1, m5[5,i], col='gold', lwd=2)

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
## This now gives the correct values

## Aah, this misses one CG; presumably the last one. Check:
data.frame(dn=kmer.seq(2), count=kmer.counts(paste0(test.monomer, "GG"), 2)$counts)
## That gives us 2 CG; but only one GG. Hence the last nucleotide is missing
## from the count. I need to fix the end condition for the function.

## test.monomer is 10 nucleotides long:
test.wc <- window.kmer.dist( test.monomer, kmer.seq(2), 6, freq=FALSE, ret.flag=1 )
## The windows are:
## CGCCAA    CG, GC, CC, CA, AA
##  GCCAAT       GC, CC, CA, AA, AT
##   CCAATG          CC, CA, AA, AT, TG
##    CAATGC             CA, AA, AT, TG, GC
##     AATGCG                AA, AT, TG, GC, CG
## test.wc$dist, is a matrix with one column per k-mer counted.
## the rows indicate the number of windows with a given count of the k-mer
## with the first row representing 0 instances and the last one being
## wsize instances. But note that the actual maximum possible is 1 + wsize - k
## I expect to have the following window counts:
##     0  1  2  3  4  5
## CG: 3  2  0  0  0  0
## GC: 1  4  0  0  0  0
## CC: 2  3  0  0  0  0
## CA: 1  4  0  0  0  0
## AA: 0  5  0  0  0  0
## AT: 1  4  0  0  0  0
## TG: 2  3  0  0  0  0
##
## The scores look correct, but, the offset is one after the beginning of the kmer

test.dimer <- paste0(rep(test.monomer, 2), collapse="")
## The windows will now be:
## CGCC AATGCG CGCCAATGCG
##                     
## The windows are:
## CGCCAATGCGCGCCAATGCG
## 
## CG, GC, CC, CA, AA,
##     GC, CC, CA, AA, AT
##         CC, CA, AA, AT, TG
##             CA, AA, AT, TG, GC
##                 AA, AT, TG, GC, CG
##                     AT, TG, GC, CG, GC, ##*##
##                         TG, GC, CG, GC, CG, ##*##
##                             GC, CG, GC, CG, GC, ##*##
##                                 CG, GC, CG, GC, CC, ##*##
##                                     GC, CG, GC, CC, CA, ##*##
##                                         CG, GC, CC, CA, AA, 
##                                             GC, CC, CA, AA, AT, 
##                                                 CC, CA, AA, AT, TG,
##                                                     CA, AA, AT, TG, GC,
##                                                         AA, AT, TG, GC, CG,

## The new windows are indicated Apart from these, the counts
## will simply double to:
## I expect to have the following window counts:
##     0  1  2  3  4  5
## CG: 6  6  3  0  0  0
## GC: 2  8  4  1  0  0
## CC: 7  8  0  0  0  0
## CA: 6  9  0  0  0  0
## AA: 5  10 0  0  0  0
## AT: 6  9  0  0  0  0
## TG: 7  8  0  0  0  0

test.2.wc <- window.kmer.dist( test.dimer, kmer.seq(2), 6, freq=FALSE, ret.flag=1 )
### That works as it should.. That suggests that the function does
### what it should.
## $dist
##      AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG
## [1,]  5 15  6 15  6  7 15  6 15 15 15  7 15  2 15 15
## [2,] 10  0  9  0  9  8  0  6  0  0  0  8  0  8  0  0
## [3,]  0  0  0  0  0  0  0  3  0  0  0  0  0  4  0  0
## [4,]  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0
## [5,]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## [6,]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## [7,]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0


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

