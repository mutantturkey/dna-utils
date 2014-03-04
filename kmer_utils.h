#include <unordered_map>
using namespace std;

// Kmer functions
void convert_kmer_to_num(char *str, const unsigned long length);
unsigned long num_to_index(const char *str, const int kmer, const long error_pos, long long *current_position);
char *index_to_kmer(unsigned long long index, long kmer);

// Utility functions
size_t strnstrip(char *s, int c, size_t len);
unsigned long long pow_four(unsigned long long x);

// Variables
typedef unordered_map<size_t,unsigned long long> kmer_map;

// file loading functions
kmer_map *get_sparse_kmer_counts_from_filename(const char *fn, const unsigned int kmer);
kmer_map *get_sparse_kmer_counts_from_file(FILE *fh, int kmer);

unsigned long long * get_dense_kmer_counts_from_file(FILE *fh, const unsigned int kmer);


void print_kmer(unsigned long long *counts, bool label, bool nonzero, unsigned int kmer);
void print_kmer(kmer_map *counts, bool label, bool nonzero, unsigned int kmer);

unsigned char alpha[256]; 
