// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "kmer_utils.h"

int main(int argc, char **argv) {
  char *filename = NULL;
  FILE *fh;

  unsigned int kmer = 0;

  bool nonzero = 0;
  bool label = 0;
  bool kmer_set = 0;

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
      filename = optarg;
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
      printf("help-text\n");
      exit(EXIT_SUCCESS);
    default:
      break;
    }
  }
  if(filename == NULL) {
    fprintf(stderr, "not input file specified with -i, reading from stdin\n");
		fh = stdin;
  }
	else {
		fh = fopen(filename, "r");
		if(fh == NULL) {
			fprintf(stderr, "Could not open %s - %s\n", filename, strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
  if(kmer == 0 && !kmer_set) {
    fprintf(stderr, "Error: kmer (-k) must be supplied\n");
    exit(EXIT_FAILURE);
  }
  if(kmer == 0) {
    fprintf(stderr, "Error: invalid kmer - '%d'.\n", kmer);
    exit(EXIT_FAILURE);
  }

  width = kmer_pow_four(kmer);

  unsigned long long *counts = kmer_counts_from_file(fh, kmer);

  // If nonzero is set, only print non zeros
  if(nonzero) {
    // if labels is set, print out our labels
    if(label) {
      for(i = 0; i < width; i++)
        if(counts[i] != 0) {
          char *kmer_str = kmer_index_to_kmer(i, kmer);
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
        if(counts[i] != 0) {
          char *kmer_str = kmer_index_to_kmer(i, kmer);
          fprintf(stdout, "%s\t%llu\n", kmer_str, counts[i]);
          free(kmer_str);
        }
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
