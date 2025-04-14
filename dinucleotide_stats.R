## These functions are somewhat specific to dinucleotide statistics
## that are not generally useful, and it is likely that better solutions
## to the problems addressed have been published. However, I have not
## found anything useful at this point in time.

## The main question that is provided here is:
##
## The expected distribution of a set of overlapping windowed
## dinucleotide counts.
##
## This seems like it should be straightforward, but it is
## more complex than one would expect because of interaction
## between overlapping di-mers. I have not addressed the
## issue of k-mers in general as the complexity of the problem
## gets worse with k. It may be possible to extend the reasoning
## used here, to k-mers, but I'm not convinced that it would
## be worth it.

## This function evaluates the graph of dimer- match and mismatches
## calculating the transition from matches to mis-matches in all
## directions. It was written in order to explore the problem domain.
##
## mnf: a named vector giving the frequencies of nucleotides
## nucs: The two nucleotides making up the dinucleotide
## depth: the number of selections to recurse to. Note that
##        the function may create up to 2^(depth+1) nodes.
##        So don't specify values of more than something like 20
##        (which is already very slow with homo-dimers).
## freq: report counts or frequencies at different nodes.
dimer.transition.graph <- function(mnf, nucs, depth, freq=FALSE){
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


## takes a graph defined as above and visualises the counts.
draw.dn.tr.graph <- function(graph, cex=1, adj.y=0.3){
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

### extracts distributions for window-lengths included in the
### graph
graph.to.distribution <- function(graph, x){
    graph <- graph[ graph$x == x, ]
    counts <- with(graph, rbind( cbind(match.n, mm), cbind(match.n+1, m) ))
    tapply( c(0:x, counts[,2]), c(rep(0,x+1), counts[,1]), sum )
}

## This calculates the expected distribution using equations
## that embody the graph. In essence it calculates the distribution
## for window of 1 dimer; then 2,..,w-1,w dimers by using the probabilities
## of transition to the specified dimer.
## The transition probability is influenced by the number of negative positions
## preceeding it (the path taken through the graph), and so we need to track
## at which points in the traversal that matches were encountered.
## This is highly un-intuitive, but the distributions obtained match distributions
## obtained from random sequences.
## the arguments to the function are as above.
windowed.dimer.exp <- function(mnf, w, k.nucs){
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
    m0.f <- 1  ## this is the initial frequency of the null branch (0 matches).
    ## we have a table of frequencies for different counts at different ages
    ## the rows indicate ages; we have an extra row and column
    counts <- matrix(0, nrow=w+1, ncol=w+1)
    counts.tmp <- counts
    ## rows indicate ages..
    tr.p[1] <- m2m.p ## this should be remove later: this is only for testing purposes
    for(i in 1:w){
        m.new <- mp.init * m0.f ## the fraction of new matches
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
