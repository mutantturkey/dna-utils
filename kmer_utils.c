// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unordered_map>

using namespace std;

typedef unordered_map<size_t,unsigned long long> kmer_map;

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

unsigned long long pow_four(unsigned long long x) {
	return (unsigned long long)1 << (x * 2);
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
unsigned long long load_specific_mers_from_file(char *fn, unsigned int kmer, size_t width, size_t *arr) { 
		FILE *fh;
		size_t arr_pos = 0;
		char line[64];

		fh = fopen(fn, "r");
		if(fh == NULL) {
			fprintf(stderr, "Error opening %s - %s\n", fn, strerror(errno));
			exit(EXIT_FAILURE);
		}	

   	while (fgets(line, sizeof(line), fh) != NULL) {
			size_t i;
			size_t len;

			len = strlen(line);
			if(len == 0)
				continue;

			len--;
			line[len] = '\0';


			if(len != kmer)  {
				fprintf(stderr, "SKIPPING: '%s' is not length %u\n", line, kmer);
				continue;
			}
			
			for(i = 0; i < len; i++) {
				line[i] = alpha[(int)line[i]];
			}
				
		 	size_t mer = num_to_index(line, kmer, width, NULL);
			if(mer == width) {
				fprintf(stderr, "SKIPPING: '%s' is a unrecognized mer\n", line);
				continue;
			}
			else {
				arr[arr_pos] = mer;
				arr_pos++;
			}
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
	if(ret == NULL || num_array == NULL) 
		exit(EXIT_FAILURE);
		

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

kmer_map *get_sparse_kmer_counts_from_file(FILE *fh, int kmer) {

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	long long i = 0;
	long long position = 0;

 	kmer_map *counts = new kmer_map();
	//
	// width is 4^kmer  
	const unsigned long long width = pow_four(kmer); 

	while ((read = getdelim(&line, &len, '>', fh)) != -1) {
		size_t k;
		char *seq;
		size_t seq_length;

		// find our first \n, this should be the end of the header
		seq = strchr(line, '\n');	
		if(seq == NULL) 
			continue;

		// point to one past that.
		seq = seq + 1;

		// strip out all other newlines to handle multiline sequences
		strnstrip(seq, '\n', strlen(seq));
		seq_length = strlen(seq);

		// relace A, C, G and T with 0, 1, 2, 3 respectively
		// everything else is 5 
		for(k = 0; k < seq_length; k++)
			seq[k] = alpha[(int)seq[k]];
		

		// loop through our seq to process each k-mer
		for(position = 0; position < (signed)(seq_length - kmer + 1); position++) {
			unsigned long long mer = 0;
			unsigned long long multiply = 1;

			// for each char in the k-mer check if it is an error char
			for(i = position + kmer - 1; i >= position; i--){
				if(seq[i] == 5) { 
					mer = width;
					position = i;
					goto next;
				}

				// multiply this char in the mer by the multiply
				// and bitshift the multiply for the next round
				mer += seq[i] * multiply;
				multiply = multiply << 2;
			}
			// use this point to get mer of our loop
			next:
			{
				// bump up the mer value in the counts array
				(*counts)[mer]++;
			}
		}
	} 

  free(line);
	fclose(fh);

	return counts;

}

unsigned long long * get_dense_kmer_counts_from_file(FILE *fh, const unsigned int kmer) {

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	long long i = 0;
	long long position = 0;

	// width is 4^kmer  
	const unsigned long long width = pow_four(kmer); 
	
	unsigned long long *counts = (unsigned long long *) calloc(width+1, sizeof(unsigned long long));
	if(counts == NULL) {
		exit(EXIT_FAILURE);
	}
	while ((read = getdelim(&line, &len, '>', fh)) != -1) {
		size_t k;
		char *seq;
		size_t seq_length;

		// find our first \n, this should be the end of the header
		seq = strchr(line, '\n');	
		if(seq == NULL) 
			continue;

		// point to one past that.
		seq = seq + 1;

		// strip out all other newlines to handle multiline sequences
		strnstrip(seq, '\n', strlen(seq));
		seq_length = strlen(seq);

		// relace A, C, G and T with 0, 1, 2, 3 respectively
		// everything else is 5 
		for(k = 0; k < seq_length; k++)
			seq[k] = alpha[(int)seq[k]];
		

		// loop through our seq to process each k-mer
		for(position = 0; position < (signed)(seq_length - kmer + 1); position++) {
			unsigned long long mer = 0;
			unsigned long long multiply = 1;

			// for each char in the k-mer check if it is an error char
			for(i = position + kmer - 1; i >= position; i--){
				if(seq[i] == 5) { 
					mer = width;
					position = i;
					goto next;
				}

				// multiply this char in the mer by the multiply
				// and bitshift the multiply for the next round
				mer += seq[i] * multiply;
				multiply = multiply << 2;
			}
			// use this point to get mer of our loop
			next:
			{
				// bump up the mer value in the counts array
				counts[mer]++;
			}
		}
	}

  //free(line);
	fclose(fh);
	return counts;

}

unsigned long long * get_dense_kmer_counts_from_filename(const char *fn, const unsigned int kmer) {
		FILE *fh = fopen(fn, "r");
		if(fh == NULL) {
			fprintf(stderr, "Could not open %s - %s\n", fn, strerror(errno));
			return 0;
		}

		return get_dense_kmer_counts_from_file(fh, kmer);
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
			for(i = 0; i < width; i++)
				if(counts->count(i) != 0) {
					char *kmer_str = index_to_kmer(i, kmer);
					fprintf(stdout, "%s\t%llu\n", kmer_str, counts->at(i));
					free(kmer_str);
				}

		}
		else {
			for(i = 0; i < width; i++)
				if(counts->count(i) != 0) 
					fprintf(stdout, "%zu\t%llu\n", i, counts->at(i));

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
