// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "kmer_utils.h"


void help() {
	printf("usage: kmer_total_count -i input_file -k kmer [-n] [-l] ...\n\n"
				 "count mers in size k from a fasta file\n"
				 "\n"
				 "  --input    -i  input fasta file to count\n"
				 "  --kmer     -k  size of mers to count\n"
				 "  --nonzero  -n  only print non-zero values\n"
				 "  --label    -l  print mer along with value\n"
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

	unsigned long long width = 0;

	unsigned long long i = 0;

	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"kmer",  required_argument, 0, 'k'},
		{"nonzero", no_argument, 0, 'n'},
		{"label", no_argument, 0, 'l'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while (1) {

		int option_index = 0;
		int c = 0;

		c = getopt_long (argc, argv, "i:k:nlvh", long_options, &option_index);

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
			case 'n':
				nonzero = true; 
				break;
			case 'l':
				label = true;
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

	width = pow_four(kmer);

	unsigned long long *counts = get_kmer_counts_from_file(fh, kmer);

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
					fprintf(stdout, "%llu\t%llu\n", i, counts[i]);

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

	free(counts);
	return EXIT_SUCCESS;
}
