// Kmer functions
void convert_kmer_to_num(char *str, const unsigned long length);
unsigned long num_to_index(const char *str, const int kmer, const long error_pos, long long *current_position);
char *index_to_kmer(unsigned long long index, long kmer);

// Utility functions
size_t strnstrip(char *s, int c, size_t len);
unsigned long long pow_four(unsigned long long x);

// Variables
const unsigned char alpha[256]; 

// file loading functions
node * get_sparse_kmer_counts_from_filename(const char *fn, const unsigned int kmer);
node * get_sparse_kmer_counts_from_file(FILE *fh, const int kmer);

unsigned long long * get_dense_kmer_counts_from_file(FILE *fh, const unsigned int kmer);


void print_dense(unsigned long long *counts, bool label, bool nonzero, unsigned int kmer);
