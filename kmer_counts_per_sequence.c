// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "kmer_utils.h"

#ifdef FORCE_GETLINE
	#include "getline/getline.h"
#endif

void help() {
	printf("usage: kmer_counts_per_sequence input_file kmer [kmer-file] ...\n\n"
				 "count mers in each sequence of size k from a fasta file\n"
				 "\n"
				 "  --input      -i  input fasta file to count\n"
				 "  --kmer       -k  size of mers to count\n"
				 "  --compliment -c  count compliment of sequences\n"
				 "  --mer-file   -m  a file containing a list of mers you are interested\n"
				 "                   in opening. this will enable output your results in\n"
				 "                   a sparse format \n"
				 "  --sparse     -s  output values in a sparse format. output is in the\n"
				 "                   order sequence_number, mer_index, value\n"
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

int main(int argc, char **argv) {


	// getdelim variables
	char *fn = NULL;
	FILE *fh = NULL;

	// for specific mers
	char *mer_fn = NULL;
	size_t num_desired_indicies = 0;
	size_t *desired_indicies = NULL;

	// for getdelim
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	// sizes
	unsigned int kmer = 0;
	size_t width = 0;

	bool sparse = false;
	bool kmer_set = false;
	bool specific_mers = false;
	bool count_compliment = false;

	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"kmer",  required_argument, 0, 'k'},
		{"compliment", required_argument, 0, 'c'},
		{"sparse", no_argument, 0, 's'},
		{"mer-file", required_argument, 0, 'm'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while (1) {

		int option_index = 0;
		int c = 0;

		c = getopt_long (argc, argv, "i:k:m:cvsh", long_options, &option_index);

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
			case 's':
				sparse = true;
				break;
			case 'm':
				mer_fn = optarg;
				specific_mers = true;
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

	if(specific_mers) {
		sparse = false;
		desired_indicies = (size_t *) malloc((width) * sizeof(size_t));
		if(desired_indicies == NULL) 
			exit(EXIT_FAILURE);
		num_desired_indicies = load_specific_mers_from_file(mer_fn, kmer, width, desired_indicies);
		if(num_desired_indicies == 0) {
			fprintf(stderr, "Error: no mers loaded from file\n"); 
			exit(EXIT_FAILURE);
		}
	}


	unsigned long long *counts = (unsigned long long *) malloc((width+ 1) * sizeof(unsigned long long));
	if(counts == NULL)  {
		fprintf(stderr, "%s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}

	unsigned long long sequence = 0;
	while ((read = getdelim(&line, &len, '>', fh)) != -1) {
		size_t k = 0;

		memset(counts, 0, width * sizeof(unsigned long long));

		// find our first \n, this should be the end of the header
		char *seq = strchr(line, '\n');	
		if(seq == NULL) 
			continue;


		// point to one past that.
		seq = seq + 1;

		// strip out all other newlines to handle multiline sequences
		size_t seq_length = strnstrip(seq, '\n', strlen(seq));

		// reset our count matrix to zero

		for(k = 0; k < seq_length; k++) {
			seq[k] = alpha[(int)seq[k]];
		}

		count_sequence(seq, seq_length, kmer, counts);

		if(count_compliment) {
			for(k = 0; k < seq_length; k++) { 
				seq[k] = compliment[(int)seq[k]];
			}
			
			reverse_string(seq, seq_length);
			count_sequence(seq, seq_length, kmer, counts);
			
		}

		if(specific_mers) {
				for(k = 0; k < num_desired_indicies; k++) {
					if(counts[desired_indicies[k]] != 0)
						fprintf(stdout, "%llu\t%zu\t%llu\n", sequence, desired_indicies[k], counts[desired_indicies[k]]);
				}
		} 
		else if(sparse) {
			for(k = 0; k < width; k++) {
				if(counts[k] != 0)
					fprintf(stdout, "%llu\t%zu\t%llu\n", sequence, k, counts[k]);
			}
		}
		else {
			for(k = 0; k < width - 1; k++)
				fprintf(stdout, "%llu\t", counts[k]);
			fprintf(stdout, "%llu\n", counts[width - 1]);
		}

		sequence++;
	}

	free(counts);
	free(line);
	free(desired_indicies);
	fclose(fh);


	return EXIT_SUCCESS;
}
