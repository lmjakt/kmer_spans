# Find spans in genomic data based on kmer-frequency data

The code in this repository explores the question as to whether we can
rapidly identify different types of genomic repeats by considering k-mer
content. The idea is rather simple:

1. Obtain a k-mer frequency spectrum.
2. Scan the genome and increment or decrement a score depending
   on the frequency of k-mers.
3. Define windows of low or high frequencies using an approach
   similar to Smith-Waterman as used in CpG-distribute.
   
Define the score at residue $i$ as:

$$
S_i = \max
	\begin{cases}
	S_{i-1} + s\\
	0
	\end{cases}
$$

where $s$ is a value that depends on the frequency of the k-mer at
residue $i$. What form the function should take is not clear, but
there are several reasonable alternatives:

$$
s = log_2(f_i / f_{med})
$$

where $f_i$ is the frequency of the k-mer at residue $i$ and $f_{med}$
is the median k-mer frequency.

We could also simply consider if the frequency is above or below some
threshold frequency $f_{t}$ (eg. the median):

$$
s = \begin{cases}
	1 & \text{if}\ f_i >= f_{t} \\
	-1 & \text{if}\ f_i < f_{t}
	\end{cases}
$$

Or use a score based on some form of a weighted rank ($r_i$) of the
kmer frequency at position $i$ and the rank at some threshold ($r_t$):

$$
s = (r_i - r_t) / r_t
$$

If we take the weighted rank of a kmer to be the proportion of the analysed sequences
containing less common k-mers, the the value will be somewhat akin to N50 values. This has
the pleasing property that the proportion of positions with positive scores will be slightly
lower than $r_t$.



# Issues and bugs

At the moment the code makes use of `qsort_r()`. This is not portable and
as I understand it, it means that it cannot be compiled on Windows or Mac
computers.
