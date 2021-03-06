// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdbool.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <getopt.h>

#include "kmer_utils.h"

template <typename array_type>
void count_sequence(const char *seq, const size_t seq_length, const unsigned int kmer, array_type *counts) {
	long long position;
	long long i;

	// loop through our seq to process each k-mer
	for(position = 0; position < (signed)(seq_length - kmer + 1); position++) {
		unsigned long long mer = 0;
		unsigned long long multiply = 1;

		// for each char in the k-mer check if it is an error char
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

void help() {
	printf("usage: kmer_total_count -i input_file -k kmer [-c] [-n] [-l] ...\n\n"
				 "count mers in size k from a fasta file\n"
				 "\n"
				 "  --input    -i  input fasta file to count\n"
				 "  --kmer     -k  size of mers to count\n"
				 "  --compliment -c  count compliment of sequences\n"
				 "  --nonzero  -n  only print non-zero values\n"
				 "  --label    -l  print mer along with value\n"
				 "  --sparse   -s  force sparse table for any mer\n"
				 "\n"
				 "Report all bugs to mutantturkey@gmail.com\n"
				 "\n"
				 "Copyright 2014 Calvin Morrison, Drexel University.\n"
				 "\n"
				 "If you are using any dna-utils tool for a publication\n"
				 "please cite your usage:\n\n"
				 "dna-utils. Drexel University, Philadelphia USA, 2014;\n"
				 "software available at www.github.com/EESI/dna-utils/\n");
}


int main(int argc, char **argv) {

	char *fn = NULL;
	FILE *fh = NULL;

	unsigned int kmer = 0;

	bool nonzero = false;
	bool label = false;
	bool kmer_set = false;
	bool count_compliment = false;
	bool force_sparse = false;

	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"kmer",  required_argument, 0, 'k'},
		{"compliment", required_argument, 0, 'c'},
		{"nonzero", no_argument, 0, 'n'},
		{"label", no_argument, 0, 'l'},
		{"sparse", no_argument, 0, 's'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while (1) {

		int option_index = 0;
		int c = 0;

		c = getopt_long (argc, argv, "i:k:cnslvh", long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
			case 'i':
				fn = optarg;
				break;
			case 'k':
				kmer = atoi(optarg);
				kmer_set = true;
				break;
			case 'c':
				count_compliment = true;
				break;
			case 'n':
				nonzero = true; 
				break;
			case 'l':
				label = true;
				break;
			case 's':
				force_sparse = true;
				break;
			case 'h':
				help();
				exit(EXIT_SUCCESS);
				break;
			case 'v':
				printf("dna-utils version " VERSION "\n");
				exit(EXIT_SUCCESS);
				break;
			default:
				break;
		}
	}
	if(argc == 1) {
		help();
		exit(EXIT_FAILURE);
	}
	if(fn == NULL) {
		fprintf(stderr, "no input file specified with -i, reading from stdin\n");
		fh = stdin;
	}
	else {
		fh = fopen(fn, "r");
		if(fh == NULL) {
			fprintf(stderr, "Could not open %s - %s\n", fn, strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
	if(!kmer_set) {
		fprintf(stderr, "Error: kmer (-k) must be supplied\n");
		exit(EXIT_FAILURE);
	}
	if(kmer == 0) { 
		fprintf(stderr, "Error: invalid kmer - '%d'.\n", kmer);
		exit(EXIT_FAILURE);
	}

	if(kmer > 12 || force_sparse) {
		kmer_map *counts = NULL;
		kmer_map *res = get_kmer_counts_from_file(counts, fh, kmer, count_compliment);

		print_kmer(res, label, nonzero, kmer);
	}
	else {
		unsigned long long *counts = NULL;
		unsigned long long *res = get_kmer_counts_from_file(counts, fh, kmer, count_compliment);
		print_kmer(res, label, nonzero, kmer);
	}

	return EXIT_SUCCESS;
}
