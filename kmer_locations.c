// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "kmer_utils.h"

void print_mer(unsigned long long mer, long const long pos, const bool labels, const bool reverse, unsigned int kmer, char *label) {
	if(labels) {
		char *kmer_str;
		bool free_ptr = false;
		if(label == NULL) {
			kmer_str = index_to_kmer(mer, kmer);
			free_ptr = true;
		}
		else {
			kmer_str = label;
		}

		if(reverse)
			printf("%s %llu c\n", kmer_str, pos);
		else
			printf("%s %llu\n", kmer_str, pos);
		if (free_ptr)
			free(kmer_str);
	}
	else {
		if(reverse)
			printf("%llu %llu c\n", mer, pos);
		else
			printf("%llu %llu\n", mer, pos);
	}
}

void print_sequence(const char *seq, // sequence
										const size_t seq_length,  // length of sequence
										const size_t global_pos,  // overall position in read
										const unsigned int kmer,  // kmer size
									  size_t *specific_mers, // specific mers  array
									  char   **specific_mers_labels, // specific mers  array
										size_t num_specific_mers,
										bool labels,  // print label instead of indicies
										bool reverse) // reverse points
{
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
		if(num_specific_mers != 0) {
			size_t j;
			for(j = 0; j < num_specific_mers; j++) {
				if(mer == specific_mers[j])
					print_mer(mer, position + global_pos, labels, reverse, kmer, specific_mers_labels[j]);
			}
		}
		else {
			print_mer(mer, position + global_pos, labels, reverse, kmer, NULL);
		}
		// skip count if error
		next: ;
	}
}

void help() {
	printf("usage: kmer_locations -i input_file -k kmer [-c] [-n] [-l] ...\n\n"
				 "print locations of mers in size k from a fasta file\n"
				 "\n"
				 "  --input      -i  input fasta file to count\n"
				 "  --kmer       -k  size of mers to count\n"
				 "  --compliment -c  count compliment of sequences (position is not in the reverse. \n"
				 "  --label      -l  print mer along with location\n"
				 "  --mer-file   -m  a file containing a list of mers you are interested\n"
				 "                   in opening. this will enable output your results in\n"
				 "                   a sparse format \n"
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

	// for getdelim
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	unsigned int kmer = 0;

	// for specific mers
	char *mer_fn = NULL;
	size_t num_desired_indicies = 0;
	size_t *desired_indicies = NULL;
	char **labeled_mers = NULL;

	bool label = false;
	bool kmer_set = false;
	bool specific_mers = false;
	bool count_compliment = false;

	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"kmer",  required_argument, 0, 'k'},
		{"compliment", required_argument, 0, 'c'},
		{"label", no_argument, 0, 'l'},
		{"mer-file", no_argument, 0, 'm'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while (1) {

		int option_index = 0;
		int c = 0;

		c = getopt_long (argc, argv, "i:k:m:clvh", long_options, &option_index);

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
			case 'm':
				specific_mers = true;
				mer_fn = optarg;
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

	if(specific_mers) {
		size_t j; 
		size_t width = pow_four(kmer);
		desired_indicies = (size_t *) malloc((width) * sizeof(size_t));
		if(desired_indicies == NULL) 
			exit(EXIT_FAILURE);

		num_desired_indicies = load_specific_mers_from_file(mer_fn, kmer, width, desired_indicies);
		if(num_desired_indicies == 0) {
			fprintf(stderr, "Error: no mers loaded from file\n"); 
			exit(EXIT_FAILURE);
		}

		labeled_mers = (char**) malloc(sizeof(char **) * num_desired_indicies);
		check_null_ptr(labeled_mers, NULL);
		for(j = 0; j < num_desired_indicies; j++) {
			labeled_mers[j] = index_to_kmer(desired_indicies[j], kmer);
		}
	}

	unsigned long long global_position = 0;
	while ((read = getdelim(&line, &len, '>', fh)) != -1) {
		size_t k = 0;

		// find our first \n, this should be the end of the header
		char *seq = strchr(line, '\n');	
		if(seq == NULL) 
			continue;


		// point to one past that.
		seq = seq + 1;

		// strip out all other newlines to handle multiline sequences
		size_t seq_length = strnstrip(seq, '\n', strlen(seq));

		for(k = 0; k < seq_length; k++) {
			seq[k] = alpha[(int)seq[k]];
		}
		
		print_sequence(seq, seq_length, global_position, kmer, desired_indicies, labeled_mers, num_desired_indicies, label, false);

		if(count_compliment) {
			for(k = 0; k < seq_length; k++) { 
				seq[k] = compliment[(int)seq[k]];
			}
			reverse_string(seq, seq_length);
			size_t rev_seq_length = seq_length;

			// chomp all errors off the beginning
			while(*seq == 5)  {
			 rev_seq_length--;
			 seq++;
			}

			print_sequence(seq, rev_seq_length, global_position, kmer, desired_indicies, labeled_mers, num_desired_indicies, label, true);
		}

		if(seq_length != 0) {
			seq_length--;
		}
		global_position += seq_length;
	}

	free(line);
	free(desired_indicies);
	free(labeled_mers);
	fclose(fh);

	return EXIT_SUCCESS;
	}
