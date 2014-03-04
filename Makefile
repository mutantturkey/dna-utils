VERSION=\"0.0.2\"
CC = g++
CFLAGS = -O3 -s -mtune=native -Wall -Wextra -DVERSION=$(VERSION) -std=c++11
DESTDIR = /usr/local/


all: kmer_utils.o libkmer.so kmer_total_count kmer_counts_per_sequence

kmer_utils.o: kmer_utils.c
	$(CC) -c kmer_utils.c -o kmer_utils.o $(CFLAGS) -fPIC
libkmer.so: kmer_utils.o
	$(CC) kmer_utils.c -o libkmer.so $(CFLAGS) -shared -fPIC
kmer_total_count: kmer_utils.o kmer_total_count.c kmer_utils.h
	$(CC) kmer_utils.o kmer_total_count.c -o kmer_total_count $(CLIBS) $(CFLAGS)
kmer_counts_per_sequence: kmer_utils.o kmer_counts_per_sequence.c kmer_utils.h
	$(CC) kmer_utils.o kmer_counts_per_sequence.c -o kmer_counts_per_sequence $(CLIBS) $(CFLAGS)

clean:
	rm -vf kmer_total_count kmer_counts_per_sequence kmer_utils.o libkmer.so

debug: CFLAGS = -ggdb -Wall -Wextra -DVERSION=$(VERSION)\"-debug\"
debug: all

install: all
	@cp -vf kmer_counts_per_sequence kmer_total_count $(DESTDIR)/bin/
	@cp -vf  libkmer.so  $(DESTDIR)/lib/

