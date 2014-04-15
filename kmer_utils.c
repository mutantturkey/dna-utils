// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <unordered_map>

using namespace std;

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

typedef unordered_map<size_t, unsigned long long, kmer_noHash_hash, kmer_eq> kmer_map;

const unsigned char alpha[256] =
{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

const char reverse_alpha[4] = { 'A', 'C', 'G', 'T' };

													// compliments
													// A  C  G  T  E  E
													// T  G  C  A  E  E
const char compliment[6] = { 3, 2, 1, 0, 5, 5};


unsigned long long pow_four(unsigned long long x) {
	return (unsigned long long)1 << (x * 2);
}

void check_null_ptr(void *ptr, const char *error) {
	if (ptr == NULL) {
		if(error != NULL)  {
			fprintf(stderr, "Error: %s - %s\n", error, strerror(errno));
		}
		else {
			fprintf(stderr, "Error: %s\n", strerror(errno));
		}
		exit(EXIT_FAILURE);
	}
}

void reverse_string(char *s, size_t len) {
	char t, *d = &(s[len - 1]);
	while (d > s) {
		t = *s;
		*s++ = *d;
		*d-- = t;
	}
}

// convert a string of k-mer size base-4 values  into a
// base-10 index
unsigned long num_to_index(const char *str, const int kmer, const long error_pos, long long *current_position) {

	int i = 0;
	unsigned long out = 0;
	unsigned long multiply = 1;

	for(i = kmer - 1; i >= 0; i--){
		if(str[i] == 5) {
			current_position += i;
			return error_pos;
		}

		out += str[i] * multiply;
		multiply = multiply << 2;
	}

	return out;
}

// return the number of loaded elements
size_t load_specific_mers_from_file(char *fn, unsigned int kmer, size_t width, size_t *arr) { 
		FILE *fh;
		size_t arr_pos = 0;
		char line[64];

		fh = fopen(fn, "r");
		check_null_ptr(fh, fn);

   	while (fgets(line, sizeof(line), fh) != NULL) {
			size_t i;
			size_t len;

			len = strlen(line);
			if(len == 0)
				continue;

			len--;
			line[len] = '\0';

			char *copy = (char *) malloc(len);
			check_null_ptr(copy, NULL);
			memcpy(copy, line, len);

			if(len != kmer)  {
				fprintf(stderr, "SKIPPING: '%s' is not length %u\n", line, kmer);
				continue;
			}
			
			for(i = 0; i < len; i++) {
				line[i] = alpha[(int)line[i]];
			}

		 	size_t mer = num_to_index(line, kmer, width, NULL);
			if(mer == width) {
				fprintf(stderr, "SKIPPING: '%s' is a unrecognized mer\n", copy);
				free(copy);
				continue;
			}
			else {
				arr[arr_pos] = mer;
				arr_pos++;
			}
			free(copy);
		}

		fclose(fh);
		return arr_pos;
}

// convert an index back into a kmer string
char *index_to_kmer(unsigned long long index, long kmer)  {

	size_t i = 0;
	size_t j = 0;
	char *num_array = (char *) calloc(kmer,  sizeof(char));
	char *ret = (char *)  calloc(kmer + 1, sizeof(char));
	check_null_ptr(num_array, NULL);
	check_null_ptr(ret, NULL);
		

	// this is the core of the conversion. modulus 4 for base 4 conversion
	while (index != 0) {
		num_array[i] = index % 4;
		index /= 4;
		i++;
	}
	
	// for our first few nmers, like AAAAA, the output would only be "A" instead
	// of AAAAA so we prepend it
	for(j = 0; j < (kmer - i); j++)
		ret[j] = 'A';

	// our offset for how many chars we prepended
	int offset = j;
	// save i so we can print it
	size_t start = i ;

	// decrement i by 1 to reverse the last i++
	i--;
	j = 0;

	// reverse the array, as j increases, decrease i
	for(j = 0; j < start; j++, i--)
		ret[j + offset] = reverse_alpha[(int)num_array[i]];
	
  // set our last character to the null termination byte
	ret[kmer + 1] = '\0';

	free(num_array);
	return ret;
}

// Strip out any character 'c' from char array 's'
// returns length of new string
size_t strnstrip(char *s, int c, size_t len) {

	size_t i = 0;
	size_t j  = 0;

	for(i = 0; i < len; i++) {
		if(s[i] != c) {
			if(j != i)
				s[j] = s[i];
			j++;
		}
	}

	s[j] = '\0';

	return j;

}

static void inc(kmer_map *map, size_t index) {
	(*map)[index]++;
}

static void inc(unsigned long long *counts, size_t index) {
	counts[index]++;
}

template <typename array_type>
void count_sequence(const char *seq, const size_t seq_length, const unsigned int kmer, array_type *counts) {
	long long position;
	long long i;

	// loop through our seq to process each k-mer
	for(position = 0; position < (signed)(seq_length - kmer + 1); position++) {
		unsigned long long mer = 0;
		unsigned long long multiply = 1;

		// for each char in the k-mer check if it is an error char

		// TODO Get rid of branch prediction - can we do it? Error flag with
		// bitshifts
		for(i = position + kmer - 1; i >= position; i--){
			if(seq[i] == 5) {
				position = i;
				goto next;
			}

			// multiply this char in the mer by the multiply
			// and bitshift the multiply for the next round
			mer += seq[i] * multiply;
			multiply = multiply << 2;
		}
		// bump up the mer value in the counts array
		inc(counts, mer);

		// skip count if error
		next: ;
	}
}


void create_array(kmer_map **counts, int kmer) {
	*counts = new kmer_map();
	(*counts)->reserve(pow_four(kmer) / 2 );
}

void create_array(unsigned long long **counts, const unsigned int kmer) {
	// width is 4^kmer  
	const unsigned long long width = pow_four(kmer); 
	
	*counts = (unsigned long long *) calloc(width+1, sizeof(unsigned long long));
	if(counts == NULL) {
		exit(EXIT_FAILURE);
	}
	else {
	}
}

template <typename array_type>
array_type * get_kmer_counts_from_file(array_type *counts, FILE *fh, const unsigned int kmer, const bool count_compliment) {

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	create_array(&counts, kmer);
	if(counts == NULL) {
		puts("Counts is null");
		exit(EXIT_FAILURE);
	}

	// TODO remove getdelim and add a fasta parsing state engine instead
	while ((read = getdelim(&line, &len, '>', fh)) != -1) {
		size_t k;
		char *seq;

		// TODO optimize out strchr with a while loop
		seq = strchr(line, '\n');	
		if(seq == NULL) 
			continue;
		// point to one past that.
		seq = seq + 1;

		// TODO can we loop jam the strnstrip and alpha? that way memory access is only done once?
		// strip out all other newlines to handle multiline sequences
		const size_t seq_length = strnstrip(seq, '\n', strlen(seq));

		// relace A, C, G and T with 0, 1, 2, 3 respectively
		// everything else is 5 
		for(k = 0; k < seq_length; k++)
			seq[k] = alpha[(int)seq[k]];
		
		count_sequence(seq, seq_length, kmer, counts);

		if(count_compliment) {
			// TODO Same for this, can we reverse and compliment at the same time?
			for(k = 0; k < seq_length; k++) { 
				seq[k] = compliment[(int)seq[k]];
			}
			
			reverse_string(seq, seq_length);
			count_sequence(seq, seq_length, kmer, counts);
			
		}
	}

	free(line);
	fclose(fh);

	return counts;
}

kmer_map * get_kmer_counts_from_filename(kmer_map *counts, const char *fn, const unsigned int kmer, const bool count_compliment) {
		FILE *fh = fopen(fn, "r");
		check_null_ptr(fh, fn);

		// Why does this work this way!!?
		kmer_map *counts2 = get_kmer_counts_from_file(counts, fh, kmer, count_compliment);
		if(counts == NULL) {
			puts("NULL IN FILENAME");
		}

		return counts2;

}


unsigned long long * get_kmer_counts_from_filename(unsigned long long *counts, const char *fn, const unsigned int kmer, const bool count_compliment) {
		FILE *fh = fopen(fn, "r");
		check_null_ptr(fh, fn);

		unsigned long long *counts2 = get_kmer_counts_from_file(counts, fh, kmer, count_compliment);
		return counts2;
}

unsigned long long * get_continuous_kmer_counts_from_file(FILE *fh, const unsigned int kmer, const bool count_compliment) {
	
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	// width is 4^kmer  
	// there's a sneaky bitshift to avoid pow dependency
	const unsigned long width = pow_four(kmer); 

	// malloc our return array
	unsigned long long * counts = (unsigned long long *) calloc((width+ 1), sizeof(unsigned long long));
	check_null_ptr(counts, NULL);

	size_t cpy_size = kmer - 1;

	char *end_of_previous = (char*) malloc(sizeof(char) * cpy_size);
	check_null_ptr(end_of_previous, NULL);
	memset(end_of_previous, 5, cpy_size);

	while ((read = getline(&line, &len, fh)) != -1) {
		if(line[0] == '>')
			continue;

		size_t j;
		size_t k;
		size_t seq_length;

		char *seq = (char *) malloc(read + cpy_size);

		memcpy(seq, end_of_previous, cpy_size);
		memcpy(seq + cpy_size, line, read);

		seq_length = read + cpy_size;

		// relace A, C, G and T with 0, 1, 2, 3 respectively
		// everything else is 5 
		for(k = cpy_size; k < seq_length; k++)
			seq[k] = alpha[(int)seq[k]];

		count_sequence(seq, seq_length, kmer, counts);

		for(j = 0, k = seq_length - (cpy_size + 1); k < seq_length; k++, j++) {
			end_of_previous[j] = seq[k];
		}

		if(count_compliment) {
			for(k = 0; k < seq_length; k++) 
				seq[k] = compliment[(int)seq[k]];

			reverse_string(seq, seq_length);
			count_sequence(seq, seq_length, kmer, counts);
		}


		free(seq);
	}

	free(end_of_previous);
	free(line);
	fclose(fh);

	return counts;
}

unsigned long long * get_continuous_kmer_counts_from_filename(const char *fn, const unsigned int kmer, const bool count_compliment) {
		FILE *fh = fopen(fn, "r");
		check_null_ptr(fh, fn);

		return get_continuous_kmer_counts_from_file(fh, kmer, count_compliment);
}

void print_kmer(unsigned long long *counts, bool label, bool nonzero, unsigned int kmer) {
	size_t width = pow_four(kmer);
	size_t i = 0;

	// If nonzero is set, only print non zeros
	if(nonzero) {
		// if labels is set, print out our labels
		if(label) {
			for(i = 0; i < width; i++)
				if(counts[i] != 0) {
					char *kmer_str = index_to_kmer(i, kmer);
					fprintf(stdout, "%s\t%llu\n", kmer_str, counts[i]);
					free(kmer_str);
				}

		}
		else {
			for(i = 0; i < width; i++)
				if(counts[i] != 0) 
					fprintf(stdout, "%zu\t%llu\n", i, counts[i]);

		}
	}
	// If we aren't printing nonzeros print everything
	else {
		if(label) {
			for(i = 0; i < width; i++) {
				char *kmer_str = index_to_kmer(i, kmer);
				fprintf(stdout, "%s\t%llu\n", kmer_str, counts[i]);
				free(kmer_str);
			} 
		}
		else {
			for(i = 0; i < width; i=i+4) {
				fprintf(stdout, "%llu\n%llu\n%llu\n%llu\n", counts[i], counts[i+1], counts[i+2], counts[i+3]);
			}
		}

	}
}

void print_kmer(kmer_map *counts, bool label, bool nonzero, unsigned int kmer) {
	size_t width = pow_four(kmer);
	size_t i = 0;

	// If nonzero is set, only print non zeros
	if(nonzero) {
		// if labels is set, print out our labels
		if(label) {
			for(auto it = counts->begin(); it != counts->end(); it++ ) {
				char *kmer_str = index_to_kmer(it->first, kmer);
				fprintf(stdout, "%s\t%llu\n", kmer_str, it->second);
				free(kmer_str);
			}
		}
		else {
			for(auto it = counts->begin(); it != counts->end(); it++ ) {
				fprintf(stdout, "%zu\t%llu\n", it->first, it->second);
			}
		}
	}
	// If we aren't printing nonzeros print everything
	else {
		if(label) {
			for(i = 0; i < width; i++) {
				char *kmer_str = index_to_kmer(i, kmer);
				if(counts->count(i) != 0)
					fprintf(stdout, "%s\t%llu\n", kmer_str, counts->at(i));
				else
					fprintf(stdout, "%s\t0\n", kmer_str);
				free(kmer_str); 
			} 
		}
		else {
			for(i = 0; i < width; i++) {
				if(counts->count(i) != 0)
					fprintf(stdout, "%llu\n", counts->at(i));
				else
					fprintf(stdout, "0\n");
			}
		}
	}
}
