dyn.load( paste(dirname(sys.frame(1)$ofile), "src/kmer_spans.so", sep="/") )

## this defines a magic number that identifies files containing
## kmer counts.
kmer.magic <- function(){ 310572L }

## this counts kmers for a single k;
## seq should be a character vector;
## the kmers will be counted for the combined
## set of sequences.
## The function returns a list containing:
## 1. nk: a vector of two elements: k, and n,
##       the total number of k-mers counted.
##       This may differ from the expected value
##       due to the presence of non-[ACTG] characters
##       in the sequence.
## 2. counts: a vector of kmer-counts
kmer.counts <- function(seq, k, with.f=TRUE){
    k <- as.integer(k)
    tmp <- .Call("kmer_counts", seq, k);
    names(tmp) <- c("n", "counts")
    tmp$n <- c(k, tmp$n)
    names(tmp$n) <- c('k', 'n')
    if(with.f)
        tmp$f <- tmp$counts / sum(tmp$counts)
    tmp
}

## identifies regions based upon arbitrary kmer scores. For example,
## this could be used to identify CpG islands by defining
## kmer.scores[CG] to be some positive integer and all other kmer-
## scores to be -1. K-mers can be weighted arbitrarily and then
## used to define regions.
## kmer.scores should be a named vector and it should be of
## length 4^k
## regions are identified using a recursively defined score:
## s[i] = max(s[i-1] + kmer_w, 0)
## where kmer_w is the weight (or score) of the current kmer.
## regions are defined from the first non-0 score to the
## maximum score as in the Smith-Waterman algorithm.
kmer.regions <- function(seq, k, kmer.scores, min.width, min.score){
    if(length(kmer.scores) != 4^k)
        stop("There should be a total of 4^k scores");
    k.seq <- kmer.seq(k)
    if(any(! (k.seq %in% names(kmer.scores))))
        stop("all kmers not defined")
    kmer.scores <- kmer.scores[ k.seq ]
    tmp <- .Call("kmer_regions_r", seq, as.integer(k), as.double(kmer.scores),
                 as.integer(min.width), as.double(min.score))
    names(tmp) <- c("n", "counts", "pos", "score")
    tmp
}

## identifies repeated spans in sequences.
## The function first calculates the k-mer spectrum of the
## sequences provided;
## it then scans the sequences, incrementing or decrementing
## a score depending on the frequency of k-mers; high frequency
## k-mers increase the score, and low frequency k-mers decrement.
## Spans with high frequency k-mers are taken from the position
## where the score starts to increase to the maximum score.
## Negative scores are not allowed, making this analogous to
## the Smith-Waterman algorithm.
## The threshold is the percentile above which scores are
## incremented. The score is incremented by:
## (kmer quantile) - thr
## where (kmer quantile) is the the proportion of kmers with 
## frequencies lower than the current k-mer.
## if thr is set to 0.5, then any k-mer more common than the
## median will increase the score; this is probably too lenient
## and the default thr value is set to 0.75 here.
kmer.low.comp.regions <- function(seq, k, min.w, min.score, thr=0.75){
    tmp <- .Call("kmer_low_comp_regions", seq, as.integer(k),
                 as.integer(min.w), as.double(min.score), thr)
    names(tmp) <- c("n", "counts", "w.rank", "pos", "score")
    tmp$pos <- t(tmp$pos)
    tmp$score <- t(tmp$score)
    tmp
}

## returns a vector of the kmers represented by the counts from
## other functions here. These are ordered by: A, C, T, G
## resulting from the binary encoding of nucleotide identities used.
kmer.seq <- function(k){
    .Call("kmer_seq_r", as.integer(k))
}

lr.regions <- function(seq, params, kmers, kmer.scores, trans.scores){
    tmp <- .Call("tr_lr_regions_r",
                 seq, as.integer(params), kmers,
                 as.double(kmer.scores), as.double(trans.scores))
    names(tmp) <- c("kmer.scores", "pos", "scores")
    rownames(tmp$kmer.scores) <- kmer.seq( as.integer(params[1]) )
    rownames(tmp$pos) <- c("seq.i", "beg", "end")
    rownames(tmp$score) <- c("score", "null")
    tmp2 <- list(kmer.scores=tmp[[1]], reg=data.frame( t(tmp$pos), t(tmp$scores) ))
    colnames(tmp2$reg)[4:5] <- c("score", "null")
    tmp2
}

## count the number of occurences of defined words of length k in sliding
## windows. Returns the distributions of those counts across the sequences
## provided.
window.kmer.dist <- function(seq, kmers, window, freq=TRUE, ret.flag=0L){
    if(length(table(nchar(kmers))) != 1)
        stop("All kmers must be of the same size")
    dists <- .Call("windowed_kmer_count_distributions_r", seq, kmers,
                   nchar(kmers[1]), as.integer(window), as.integer(ret.flag))
    names(dists) <- c("dist", "seq.i", "scores")
    colnames(dists$dist) <- kmers
    if(!is.null(dists$scores)){
        for(i in 1:length(dists$scores)){
            colnames(dists$scores[[i]]) <- kmers
        }
    }
    if(freq)
        dists$dist <- dists$dist / colSums(dists$dist)
    dists
}

## DEPENDANCY WARNING: This function currently uses Biostrings to read in sequences
## it will be rewritten in the future to remove this dependency and to allow for
## more efficient throughput.
## this will simply read in a load of sequence files and then output counts to a file
## specified by the prefix.
## The files will be binary encoded as 32 bit integers.
## each file can hold data from several k-mers; The file will contain (offsets given in bytes)
## 0: signed integer: magic number identifying the file as correct
## 4: signed integer: number of k's (n) counted
## 8: signed integer: number of k-mers for each k
## 8 + n * 4: signed integers: the counts for the first k-mer.
##
## the actual values of k can be derived as: log2(length) / 2.
## returns a list of in and output files and information about the sequence
kmers.to.file <- function(seq.f, out.prefix, k, min.l=1e5, magic=kmer.magic()){
    out.f <- paste0(out.prefix, "counts_", paste(k, collapse="_"), ".bin")
    seq.size <- 0
    seq.fsize <- 0
    seq.fl <- 0
    read.count <- function(){
        seq <- readDNAStringSet(seq.f)
        seq.size <<- sum(nchar(seq))
        seq <- seq[ nchar(seq) >= min.l ]
        seq.fsize <<- sum(nchar(seq))
        seq.fl <<- length(seq)
        if(length(seq) < 1)
            stop("No sequence after length filtering")
        seq <- as.character(seq)
        counts <- lapply(k, function(x){
            kmer.counts(seq, x)
        })
        counts
    }
    counts <- try(read.count(), silent=TRUE )
    if(class(counts) == "try-error")
        return(list(seq.f, NA, seq.size, seq.fsize, seq.fl))
    con <- file(out.f, open="wb")
    writeBin(as.integer(magic), con)
    writeBin(as.integer(length(counts)), con )
    for(x in counts)
        writeBin(as.integer(length(x$counts)), con)
    for(x in counts)
        writeBin(as.integer(x$counts), con)
    close(con)
    list(seq.f, out.f, seq.size, seq.fsize, seq.fl)
}

read.kmers <- function(fname, magic=kmer.magic()){
    con <- file(fname, open="rb")
    m1 <- readBin( con, "integer", n=1 )
    if(m1 != magic){
        close(con)
        return(FALSE)
    }
    kn <- readBin(con, "integer", n=1)
    if(kn < 1){
        close(con)
        return(FALSE)
    }
    ks <- readBin(con, "integer", n=kn)
    counts <- lapply(ks, function(n){
        readBin(con, "integer", n=n)
    })
    close(con)
    list(k=as.integer(log2(ks) / 2), counts=counts)
}
