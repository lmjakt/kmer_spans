#define _GNU_SOURCE
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for qsort

/*
  kmer counts will be kept in a single array with 
  the kmer index simply being the two-bit representation
  of A, C, G, T obtained from concatenation of bits 2 and 3
  of the individual bytes. If c is the nucleotide and b is
  the two bit representation:
  
  b <- (c >> 1) & 3
  This maps as:
  A: 0, C: 1, G: 3, T: 2.
  
  Hence to update the kmer offset represented as an unsigned
  integer (or possibly long) we can simply do
  
  offset <- (offset << 2) | ((c >> 1) & 3)
  
  and to update the counts:
  
  counts[offset]++

  Note that we will check for Ns in the sequence and not determine an
  offset for any words that contains an N;
*/

// a macro to update an offset as defined here
#define UPDATE_OFFSET(off, c) ( ((off) << 2) | (( (c) >> 1) & 3) )
#define LC(c) ( (c) | 0x20 )
#define MAX_K 16
#define REG_N 100

// Used to translate offsets to sequences.
const char NUC[4] = {'A', 'C', 'T', 'G'};

// Holds information about n regions
// with space for up to capacity regions
// Note that int are stored as a single block
struct seq_regions {
  // data is stored in contiguous blocks
  // one for doubles and one for ints
  // with the data mixed as in R.
  // int_data stores seq_id, start, end
  int *int_data;
  size_t ints_nrow;
  // double_data score and entropy
  size_t doubles_nrow;
  double *double_data;
  size_t n;
  size_t capacity;
};

void seq_regions_free(struct seq_regions *sr){
  free(sr->int_data);
  free(sr->double_data);
  sr->capacity = 0;
  sr->n = 0;
}

struct seq_regions init_regions(size_t i_size){
  if(i_size <= 0)
    i_size = REG_N;
  struct seq_regions sr;
  sr.capacity=i_size;
  sr.ints_nrow = 3;
  sr.doubles_nrow = 2;
  sr.n = 0;
  sr.int_data = calloc(sizeof(int), sr.capacity * sr.ints_nrow);
  sr.double_data = calloc(sizeof(double), sr.capacity * sr.doubles_nrow);
  return(sr);
}

void seq_regions_push(struct seq_regions *sr,
		      int seq_id, int start, int end,
		      double score, double entropy){
  size_t n = sr->n;
  if(n == sr->capacity){
    sr->capacity *= 2;
    int *old_ints = sr->int_data;
    double *old_doubles = sr->double_data;
    sr->int_data = calloc(sizeof(int), sr->capacity * sr->ints_nrow);
    sr->double_data = calloc(sizeof(double), sr->capacity * sr->doubles_nrow);
    memcpy(sr->int_data, old_ints, n * sizeof(int) * sr->ints_nrow);
    memcpy(sr->double_data, old_doubles, n * sizeof(double) * sr->doubles_nrow);
    free( old_ints );
    free( old_doubles );
  }
  sr->int_data[ n * sr->ints_nrow ] = seq_id;
  sr->int_data[ 1 + n * sr->ints_nrow ] = start;
  sr->int_data[ 2 + n * sr->ints_nrow ] = end;
  sr->double_data[ n * sr->doubles_nrow ] = score;
  sr->double_data[ 1 + n * sr->doubles_nrow ] = entropy;
  sr->n++;
}

// allocates memory for a counts table;
// Not needed if allocated using R functions.
int* init_counts(int k){
  if(k > MAX_K)
    return(0);
  return( calloc( sizeof(int), 1 << (2*k) ));
}

// sets at_end to 1 if at end
size_t skip_n(const char *seq, size_t i){
  while(seq[i] && LC(seq[i]) == 'n')
    ++i;
  return(i);
}

// returns the value of j; if it runs out of sequence (i.e. last k-mer)
// or it finds an N. The caller must check the return value.
size_t init_kmer(const char *seq, size_t i, unsigned long *offset, int k){
  size_t j = 0;
  while(seq[i]){
    *offset = 0;
    for(j=0; j < k && seq[i+j] && LC(seq[i+j]) != 'n'; ++j){
      *offset = UPDATE_OFFSET(*offset, seq[i+j]);
    }
    if(seq[i+j] == 0 || j == k)
      break;
    // otherwise we hit Ns again;
    i = skip_n(seq, i + j);
  }
  return( i + j );
}

// returns the number of words counted.
size_t sequence_kmer_count(const char *seq, int *counts, int k){
  size_t i = 0;
  unsigned long offset = 0;
  unsigned long word_count = 0;
  size_t mask = (1 << (2*k)) - 1;
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer(seq, i, &offset, k);
    while(seq[i] && LC(seq[i]) != 'n'){
      counts[ offset & mask ]++;
      ++word_count;
      // note that UPDATE_OFFSET can handle the NULL terminator
      // so we do not need to check for it.
      offset = UPDATE_OFFSET(offset, seq[i]);
      ++i;
    }
  }
  return(word_count);
}

// sets the seq to the sequence represented by the given
// offset. Assumes that *seq is a char array of size
// k + 1. The function will set seq[k] to 0, though it
// would be more efficient to do so elsewhere.
int kmer_seq(unsigned char *seq, unsigned int k, size_t offset){
  if(k > MAX_K)
    return(-1);
  size_t mask = 3;
  seq[k] = 0;
  for(ssize_t off=k-1; off >= 0; off--){
    seq[ off ] = NUC[ (offset & mask) ];
    offset >>= 2;
  }
  return(0);
}

// a helper function for qsort. This allows the sorting of indices
// and should be used with qsort_r( )
// a and b should be size_t (hence unsigned)
// Note that we could use ksort to make this faster;
int comp_index_int(const void *a, const void *b, void *c){
  // this only works if the values do not overflow;
  // it seems exceedingly unlikely for this to happen;
  size_t i=*((size_t*) a);
  size_t j=*((size_t*) b);
  int* v = (int*) c;
  return( v[i] - v[j] );
}

// This creates a weighted rank; that is, the proportion of
// kmers in the sequences that have frequencies below the
// given kmer; for the most rare k-mer this will be 0;
void rank_kmers_w(int *counts, size_t k, double *ranks, double total_kmer_count){
  // sort indices and then assign ranks from the indices:
  size_t n_elt = (1 << (k * 2));
  size_t *indices = malloc(sizeof(size_t) * n_elt);
  //  size_t *ranks = malloc(sizeof(size_t) * n_elt);
  for(size_t i=0; i < n_elt; ++i)
    indices[i] = i;
  qsort_r( (void*)indices, n_elt, sizeof(size_t), 
	   &comp_index_int, counts );
  ranks[0] = 0;
  for(size_t i=1; i < n_elt; ++i)
    ranks[ indices[i] ] = ranks[indices[i-1]] + (counts[indices[i-1]] / total_kmer_count);
  free(indices);
}

// it would be better if we specify a scoring function here; however,
// we will define that at a later point. For now we try something out
// in order to define the structure of the program.
// here we will use ranks; these should be supplied by the caller
// as they are common to the complete analysis.
void low_complexity_regions(const char *seq, int seq_id, size_t start, size_t min_width,
			    double min_score, double *kmer_ranks, int k, double threshold,
			    struct seq_regions *reg){
  // first define the median of the kmer_spectrum;
  //   size_t ks = (1 << (2*k)); (not used)
  // kmer_ranks will be returned
  //  size_t *kmer_ranks = rank_counts(kmer_spectrum, k);
  // double med_rank = 0.5; // (double)(ks / 2) - (ks % 2 == 0 ? 0.5 : 0);
  // and then start calculating offsets.
  // This is very similar to when we count kmers, but..
  unsigned long offset = 0;
  size_t mask = (1 << (2*k)) - 1;
  double score = 0;
  double last_score = 0;
  double max_score = 0;
  size_t reg_beg = 0;
  size_t max_pos = 0;
  size_t i = start;
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer(seq, i, &offset, k);
    last_score = max_score = score = 0;
    while(seq[i] && LC(seq[i]) != 'n'){
      double s = (kmer_ranks[ mask & offset ] - threshold); // / threshold;
      score = last_score + s;
      score = score > 0 ? score : 0;
      if(last_score == 0 && score > 0){
	reg_beg = i;
	max_pos = i;
	max_score = score;
      }
      if(score == 0 && last_score > 0){
	// then we should consider adding a region to regions
	// add_region..
	if(max_pos - reg_beg >= min_width && max_score >= min_score){
	  seq_regions_push(reg, seq_id, reg_beg, max_pos, max_score, 0);
	  i = 1 + max_pos -k;
	  break;
	}
	max_score = 0;
	max_pos = i;
      }
      if(score > max_score){
	max_score = score;
	max_pos = i;
      }
      last_score = score;
      offset = UPDATE_OFFSET(offset, seq[i]);
      ++i;
      // note that UPDATE_OFFSET can handle the NULL terminator
      // so we do not need to check for it.
    }
    // check for end region or last region before Ns
    if(score > 0){
      // then we should consider adding a region to regions
      // add_region..
      if(max_pos - reg_beg >= min_width && max_score >= min_score){
	seq_regions_push(reg, seq_id, reg_beg, max_pos, max_score, 0);
	i = 1 + max_pos - k;
      }
    }
  }
}

// seq: a sequence to be scanned
// seq_id: an integer identifying the sequence
// k: the kmer used
// kmer_scores:  scores associated with kmers;
//               these are used to define the starting score for
//               any blocks.
// trans_scores:  scores for transitions that give the
//                probability of next character given
//                the preceding (k-1) characters.
//              The rows of these should be arranged according
//              to the 2-bit representation used here, such that
//              the terminal two bits (LSB) give the emitted nucleotide.
//              Each will be assumed to have k^4 rows, but this will not
//              be checked. (Nucleotides will be recycled as, A, C, T, G)
// regions : a pointer to a seq_regions block. The regions identified will
//           be added to this block. For this, we need the seq_id.
//           Note that seq_regions, stores the seq_id, begin and end
//           positions as well as two doubles. The two doubles have
//           been used to store score and entropy. We will give the score for one.
void find_kmer_tr_lr_regions(const char *seq, int seq_id,
			     int k,
			     const double *kmer_scores,
			     const double *trans_scores,
			     struct seq_regions *regions,
			     int min_region_length){
  size_t i = 0;
  unsigned long offset = 0;
  size_t mask = (1 << (2*k)) - 1;
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer(seq, i, &offset, k);
    if(!seq[i] || !seq[i+1])
      break;
    offset = offset & mask;
    double last_score = 0;
    double max_score = 0;
    int max_score_pos = 0;
    int reg_begin = 0;
    double score = kmer_scores[offset];
    //    Rprintf("init: %c\t%ld\t%ld : %.3f\n", seq[i], i, offset, score);
    score = score < 0 ? 0 : score;
    // unfortunately I do not see any way around repeating this bit of code
    // I should try to refactor it at some later point.
    if( score > 0){
      max_score = score;
      max_score_pos = i;
      reg_begin = i;
    }
    last_score = score;
    //++i; the initialisation takes care of this
    while(seq[i] && LC(seq[i]) != 'n'){
      offset = mask & UPDATE_OFFSET(offset, seq[i]);
      score = last_score + trans_scores[offset];
      if(score > max_score){
	max_score = score;
	max_score_pos = i;
      }
      //      Rprintf("trans: %c\t%ld\t%ld : %.3f  max: %.3f\n", seq[i], i, offset, score, max_score);
      score = score < 0 ? 0 : score;
      if(last_score == 0 & score > 0){
	max_score = score;
	max_score_pos = i;
	reg_begin = i;
      }
      if(score == 0 && last_score > 0){
	if(max_score_pos - reg_begin >= min_region_length){
	  // add a region to the regions.. and since this is for R I will
	  // add one to the positions.
	  seq_regions_push(regions, seq_id, 1+reg_begin, 1+max_score_pos, max_score, 0);
	}
	i = max_score_pos;
	score = last_score = max_score = 0;
	reg_begin = i;
	max_score_pos = 0;
	// we also have to reinitialise the k-mer as we are jumping back in the sequence here:
	init_kmer(seq, 1+i-k, &offset, k);
	//continue; // we don't want to increment i as this is done by the init_kmer function.
      }
      last_score = score;
      ++i;
    }
    // check for a terminal region:
    if(max_score > 0 && max_score_pos - reg_begin >= min_region_length)
      seq_regions_push(regions, seq_id, 1+reg_begin, 1+max_score_pos, max_score, 0);
  }
}

// Return frequency spectrum to an R-session
SEXP kmer_counts(SEXP seq_r, SEXP k_r){
  if(TYPEOF( seq_r ) != STRSXP || length(seq_r) < 1 )
    error("seq_r must be a character vector of length at least one");
  if(TYPEOF( k_r ) != INTSXP || length(k_r) < 1 )
    error("k_r must be an integer vector of length at least one");
  int seq_n = length(seq_r);
  //  int k_n = length(k_r);
  int k = INTEGER( k_r )[0];
  if(k < 1 || k > MAX_K)
    error("k must be a positive integer less than 1+MAX_K");

  // return the vector of counts and the total number of
  // words counted (as a double)
  size_t counts_size = 1 << (2 * k);
  SEXP ret_value = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT( ret_value, 0, allocVector(REALSXP, 1));
  SET_VECTOR_ELT( ret_value, 1, allocVector(INTSXP, counts_size));
  double *n_counts = REAL(VECTOR_ELT(ret_value, 0));
  *n_counts = 0;
  int *counts = INTEGER(VECTOR_ELT(ret_value, 1));
  memset(counts, 0, sizeof(int) * counts_size);

  for(int i=0; i < seq_n; ++i){
      SEXP seq = STRING_ELT(seq_r, i);
      int seq_l = length(seq);
      if(seq_l < k)
	continue;
      size_t n = sequence_kmer_count(CHAR(seq),
				     counts,
				     k);
      n_counts[0] = n_counts[0] + (double)n;
  }
  UNPROTECT(1);
  return(ret_value);
}

SEXP kmer_low_comp_regions(SEXP seq_r, SEXP k_r, SEXP min_width_r, SEXP min_score_r, SEXP threshold_r){
  if(TYPEOF( seq_r ) != STRSXP || length(seq_r) < 1 )
    error("seq_r must be a character vector of length at least one");
  if(TYPEOF( k_r ) != INTSXP || length(k_r) < 1 )
    error("k_r must be an integer vector of length at least one");
  if(TYPEOF( min_width_r ) != INTSXP || length(min_width_r) != 1 )
    error("the minimum width must be an integer vector of length 1");
  if(TYPEOF( min_score_r ) != REALSXP || length(min_score_r) != 1 )
    error("the minimum score must be a REAL vector of length 1");
  if(TYPEOF( threshold_r ) != REALSXP || length(threshold_r) != 1 )
    error("the threshold must be a REAL vector of length 1");

  
  int seq_n = length(seq_r);
  //  int k_n = length(k_r);
  int k = INTEGER( k_r )[0];
  int min_width = INTEGER(min_width_r)[0];
  double min_score = REAL(min_score_r)[0];
  double threshold = REAL(threshold_r)[0];
  if(threshold <= 0 || threshold >= 1)
    error("the threshold must be between 0 and 1");

  // First obtain the spectrum for the full set of sequences;
  // This repeats code in the kmer_counts function. This is unfortunate, 
  // and we should refactor this code.
  size_t counts_size = 1 << (2 * k);
  SEXP ret_value = PROTECT(allocVector(VECSXP, 5));
  // the elements of ret value are:
  // 1. a vector containing the number of kmers counted and the number of positive scores
  // 2. the kmer counts
  // the k-mer ranks and the region integers and region doubles
  SET_VECTOR_ELT( ret_value, 0, allocVector(REALSXP, 2));
  double *n_counts = REAL(VECTOR_ELT(ret_value, 0));
  *n_counts = 0;
  SET_VECTOR_ELT( ret_value, 1, allocVector(INTSXP, counts_size));
  int *counts = INTEGER(VECTOR_ELT(ret_value, 1));
  memset(counts, 0, sizeof(int) * counts_size);
  // note that we will have overflows if more than 2^31 (i.e. k=15)
  SET_VECTOR_ELT( ret_value, 2, allocVector(REALSXP, counts_size));
  double *kmer_ranks = REAL(VECTOR_ELT(ret_value, 2));
  // pos score counts the number of k-mers that have a positive score
  // and reports this back. This is to ensure that the statistical assumptions
  // are correct;

  for(int i=0; i < seq_n; ++i){
    SEXP seq = STRING_ELT(seq_r, i);
    int seq_l = length(seq);
    if(seq_l < k)
      continue;
    size_t n = sequence_kmer_count(CHAR(seq),
				   counts,
				   k);
    n_counts[0] = n_counts[0] + (double)n;
  }
  rank_kmers_w(counts, k, kmer_ranks, *n_counts);
  // REG_N should probably be an option.
  struct seq_regions regions = init_regions(REG_N);
  for(int i=0; i < seq_n; ++i){
    SEXP seq = STRING_ELT(seq_r, i);
    int seq_l = length(seq);
    if(seq_l < k)
      continue;
    low_complexity_regions(CHAR(seq), i, 0, min_width, min_score,
			   kmer_ranks, k, threshold, &regions);
  }
  n_counts[1] = 0;
  SET_VECTOR_ELT(ret_value, 3, allocMatrix(INTSXP, regions.ints_nrow, regions.n));
  SET_VECTOR_ELT(ret_value, 4, allocMatrix(REALSXP, regions.doubles_nrow, regions.n));
  memcpy( INTEGER(VECTOR_ELT(ret_value, 3)), regions.int_data, sizeof(int) * regions.ints_nrow * regions.n);
  memcpy( REAL(VECTOR_ELT(ret_value, 4)), regions.double_data, sizeof(double) * regions.doubles_nrow * regions.n);
  seq_regions_free(&regions);
  UNPROTECT(1);
  return(ret_value);
}

SEXP kmer_seq_r(SEXP k_r){
  if(TYPEOF(k_r) != INTSXP || length(k_r) != 1)
    error("k_r should be an integer of length 1");
  unsigned int k = (unsigned int)(INTEGER(k_r)[0]);
  if(k > MAX_K || k < 1)
    error("k_r (%d) should be smaller than MAX_K (%d) and larger than 0", k, MAX_K);
  size_t n = 1 << (2 * k);
  SEXP ret_value = PROTECT(allocVector(STRSXP, n));
  unsigned char *buffer = malloc(sizeof(unsigned char) * (k+1));
  for(size_t i=0; i < n; ++i){
    kmer_seq(buffer, k, i);
    SET_STRING_ELT(ret_value, i, mkChar((const char*)buffer));
  }
  free(buffer);
  UNPROTECT(1);
  return(ret_value);
}

// seq_f: the sequences (a character vector)
// params_r: paramters (integer vector: k, min_length)
// kmers_r: the kmer sequences as a character vector
//          this is so that I can make use of kmer spectra
//          created using biostring functions; these have a different order
//          kmers_r should be of the same length as freq_a, freq_b;
// kmer_scores_r: The scores associated with beginning a region with a given kmer
// trans_scores_r: The transition scores associated extending a region by one nucleotide
SEXP tr_lr_regions_r(SEXP seq_r, SEXP params_r,
		     SEXP kmers_r,
		     SEXP kmer_scores_r, SEXP trans_scores_r){
  if(TYPEOF(seq_r) != STRSXP || length(seq_r) < 1)
    error("seq_r should be a character vector of of positive length");
  if(TYPEOF(params_r) != INTSXP || length(params_r) != 2)
    error("params_r should have two integers (k, and min_length)");
  if(TYPEOF(kmers_r) != STRSXP)
    error("kmers_r should be a character vector");
  if(TYPEOF(kmer_scores_r) != REALSXP || TYPEOF(kmer_scores_r) != REALSXP)
    error("freq_a and freq_b should be double vectors");
  int seq_n = length(seq_r);
  int k = INTEGER(params_r)[0];
  int min_length = INTEGER(params_r)[1];
  if(k < 1 || k > MAX_K)
    error("k should be a positive value less than MAX_K");
  if(min_length < 0)
    error("min_length should be a positive integer");
  size_t kmers_size = 1 << (2 * k);
  if((size_t)length(kmers_r) != kmers_size ||
     (size_t)length(kmer_scores_r) != kmers_size ||
     (size_t)length(trans_scores_r) != kmers_size)
    error("kmers_r, freq_a, freq_b should all be 4^k long");
  // we need to rearrange the order of the kmer scores in order to
  // fit with 2 bit encoding
  double *kmer_scores_in = REAL(kmer_scores_r);
  double *trans_scores_in = REAL(trans_scores_r);
  // We can return the rearranged kmer counts as part of a list
  // This should be optional though. But for now we will include
  // everything for troubleshooting.
  SEXP ret_value = PROTECT(allocVector(VECSXP, 3));
  // the four columns are, rearranged kmer spectra and the transition matrices;
  SET_VECTOR_ELT( ret_value, 0, allocMatrix(REALSXP, kmers_size, 2));
  double *spectra = REAL(VECTOR_ELT(ret_value, 0));
  double *kmer_scores = spectra;
  double *trans_scores = spectra + 1 * kmers_size;
  // We then need to set these values appropriately;
  size_t offset = 0;
  for(size_t i=0; i < kmers_size; ++i){
    int ii = init_kmer(CHAR(STRING_ELT(kmers_r, i)), 0, &offset, k);
    // ii should be checked to make sure that everything when well. It should be equal to k
    if(ii != k)
      Rprintf("ii  %d != %d\n", ii, k);
    kmer_scores[offset] = kmer_scores_in[i];
    trans_scores[offset] = trans_scores_in[i];
  }
  // The second element should contain the regions. But we do not know how many
  struct seq_regions regions = init_regions(REG_N); // 0 uses the default size
  // and go through the sequences:
  for(int i=0; i < seq_n; ++i){
    find_kmer_tr_lr_regions( CHAR(STRING_ELT(seq_r, i)), i+1,
			     k, kmer_scores, trans_scores,
			     &regions,
			     min_length);
  }
  // Here copy the regions to appropriate R structures:
  SET_VECTOR_ELT(ret_value, 1, allocMatrix(INTSXP, regions.ints_nrow, regions.n));
  memcpy( INTEGER(VECTOR_ELT(ret_value, 1)), regions.int_data, sizeof(int) * regions.ints_nrow * regions.n );
  SET_VECTOR_ELT(ret_value, 2, allocMatrix(REALSXP, regions.doubles_nrow, regions.n));
  memcpy( REAL(VECTOR_ELT(ret_value, 2)), regions.double_data, sizeof(double) * regions.doubles_nrow * regions.n );
  // 
  seq_regions_free(&regions);
  UNPROTECT(1);
  return(ret_value);
}

static const R_CallMethodDef callMethods[] = {
  {"kmer_counts", (DL_FUNC)&kmer_counts, 2},
  {"kmer_low_comp_regions", (DL_FUNC)&kmer_low_comp_regions, 5},
  {"kmer_seq_r", (DL_FUNC)&kmer_seq_r, 1},
  {"tr_lr_regions_r", (DL_FUNC)&tr_lr_regions_r, 5},
  {NULL, NULL, 0}
};

void R_init_kmer_spans(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
