VERSION=\"0.0.5\"
CC = g++
CFLAGS = -O3 -s -march=native -Wall -Wextra -DVERSION=$(VERSION) -std=c++11
DESTDIR = /usr/local/

all: libkmer.so kmer_total_count kmer_counts_per_sequence kmer_utils.o kmer_locations

kmer_utils.o: kmer_utils.c
	$(CC) -c kmer_utils.c -O $(CFLAGS) -fPIC

libkmer.so: kmer_utils.o
	$(CC) kmer_utils.o -o libkmer.so $(CFLAGS) -shared -fPIC

kmer_total_count: kmer_utils.o kmer_total_count.c kmer_utils.h
	$(CC) kmer_utils.o kmer_total_count.c -o kmer_total_count $(CLIBS) $(CFLAGS)

kmer_counts_per_sequence: kmer_utils.o kmer_counts_per_sequence.c kmer_utils.h
	$(CC) kmer_utils.o kmer_counts_per_sequence.c -o kmer_counts_per_sequence $(CLIBS) $(CFLAGS)

kmer_continuous_count: kmer_utils.o kmer_continuous_count.c kmer_utils.h
	$(CC) kmer_utils.o kmer_continuous_count.c -o kmer_continuous_count $(CLIBS) $(CFLAGS)

kmer_locations: kmer_utils.o kmer_locations.c kmer_utils.h
	$(CC) kmer_utils.o kmer_locations.c -o kmer_locations $(CLIBS) $(CFLAGS)
clean:
	rm -vf kmer_total_count kmer_counts_per_sequence kmer_continuous_count kmer_utils.o libkmer.so
	rm -vf kmer_total_count kmer_counts_per_sequence kmer_continuous_count libkmer.so kmer_utils.o

debug: CFLAGS = -ggdb -Wall -Wextra -DVERSION=$(VERSION)\"-debug\" -std=c++11
debug: all

install: all
	@cp -vf kmer_counts_per_sequence kmer_total_count kmer_continuous_count $(DESTDIR)/bin/
	@cp -vf  libkmer.so  $(DESTDIR)/lib/

