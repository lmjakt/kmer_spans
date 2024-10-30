dyn.load( paste(dirname(sys.frame(1)$ofile), "src/kmer_spans.so", sep="/") )

## this defines a magic number that identifies files containing
## kmer counts.
kmer.magic <- function(){ 720531L }

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
kmer.counts <- function(seq, k){
    k <- as.integer(k)
    tmp <- .Call("kmer_counts", seq, k);
    names(tmp) <- c("n", "counts")
    tmp$n <- c(k, tmp$n)
    names(tmp$n) <- c('k', 'n')
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
    if(m1 != magic)
        return(FALSE)
    kn <- readBin(con, "integer", n=1)
    if(kn < 1)
        return(FALSE)
    ks <- readBin(con, "integer", n=kn)
    counts <- lapply(ks, function(n){
        readBin(con, "integer", n=n)
    })
    list(k=as.integer(log2(ks) / 2), counts=counts)
}
