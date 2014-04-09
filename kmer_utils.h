// dna-util's function library.
#include <unordered_map>
using namespace std;

// Kmer functions
void convert_kmer_to_num(char *str, const unsigned long length);
unsigned long num_to_index(const char *str, const int kmer, const long error_pos, long long *current_position);
char *index_to_kmer(unsigned long long index, long kmer);

// Utility functions


// strip char 'c' out of char array *s of length len
size_t strnstrip(char *s, int c, size_t len);

// reverse char arry *s of length len
void reverse_string(char *s, size_t len);

// quicky calculate 4^x 
unsigned long long pow_four(unsigned long long x);

// check if pointer is null. a helper for dealing with NULL 
// return values as errors. Calls strerror and quits if 
// ptr is null, optionally takes *error char array as 
// a error to output
void check_null_ptr(void *ptr, const char *error);

template <typename array_type>
void count_sequence(const char *seq, const size_t seq_length, const unsigned int kmer, array_type *counts);

// Variables
typedef struct {
	size_t operator() (const size_t &k) const {
	return k;
	}
} kmer_noHash_hash;

typedef struct {
	bool operator() (const size_t &x, const size_t &y) const {
		return x == y;
	}
} kmer_eq; 

typedef unordered_map<size_t,unsigned long long, kmer_noHash_hash, kmer_eq> kmer_map;

unsigned char alpha[256]; 
unsigned char reverse_alpha[4];
unsigned char compliment[5];

// open file from filename in char array *fn, and try and parse in one mer per
// line, of size kmer, and store the indicies of those mers in the *arr
// pointer;
unsigned long long load_specific_mers_from_file(const char *fn, unsigned int kmer, size_t width, size_t *arr);

unsigned long long * get_continuous_kmer_counts_from_filename(const char *fn, const unsigned int kmer, const bool count_compliment);
unsigned long long * get_continuous_kmer_counts_from_file(FILE *fh, const unsigned int kmer, const bool count_compliment);


template <typename array_type>
array_type * get_kmer_counts_from_file(array_type *counts, FILE *fh, const unsigned int kmer, const bool count_compliment);

kmer_map           *get_kmer_counts_from_filename(kmer_map           *counts, const char *fn, const unsigned int kmer, const bool count_compliment);
unsigned long long *get_kmer_counts_from_filename(unsigned long long *counts, const char *fn, const unsigned int kmer, const bool count_compliment);


size_t load_specific_mers_from_file(char *fn, unsigned int kmer, size_t width, size_t *arr);

// print functions
void print_kmer(unsigned long long *counts, bool label, bool nonzero, unsigned int kmer);
void print_kmer(kmer_map *counts, bool label, bool nonzero, unsigned int kmer);
