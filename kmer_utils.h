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
unsigned long long load_specific_mers_from_file(const char *fn, unsigned int kmer, size_t width, size_t *arr);
node * get_kmer_counts_from_filename(const char *fn, const unsigned int kmer);
node * get_kmer_counts_from_file(FILE *fh, const int kmer);
