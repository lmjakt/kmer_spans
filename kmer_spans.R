dyn.load( paste(dirname(sys.frame(1)$ofile), "src/kmer_spans.so", sep="/") )

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
kmer.low.comp.regions(seq, k, min.w, min.score, thr=0.75){
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
