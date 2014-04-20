# dna-utils 


## building

Building is easy using make.

Linux:

    make
		make install

OS/X or if you don't have getline:

    make -f Makefile.osx
		make install
This repository contains general utilities for processing sequences in fasta files.


## Tools included ##

All of our tools use an alphebetical indexing scheme like this:

    AAAA = 0
    AAAC = 1
    AAAG = 2
    AAAT = 3
    AACA = 4
    ...
    

### kmer_total_count 

this program will count each kmer in a fasta file, and print to standard out.

#### Usage

    usage: kmer_total_count -i input_file -k kmer [-n] [-l] ...

    count mers in size k from a fasta file

      --input    -i  input fasta file to count
      --kmer     -k  size of mers to count
      --nonzero  -n  only print non-zero values
      --label    -l  print mer along with value

    Report all bugs to mutantturkey@gmail.com

    Copyright 2014 Calvin Morrison, Drexel University.

    If you are using any dna-utils tool for a publication
    please cite your usage:

    dna-utils. Drexel University, Philadelphia USA, 2014;
    software available at www.github.com/EESI/dna-utils/

#### Examples

a basic example, where we specify the k-mer size and input file.

    calvin@barnabas:~/dna-utils$ ./kmer_total_count -i SuperManSequences.fasta -k 8 
    2946
    1161
    14141
    ...

it also supports input from stdin, which is great for combining with compression programs

    calvin@barnabase:~/dna-utils$ gzip -dc ~/super_big_fasta_file.fasta.gz | ./kmer_total_count -k 8
    234523
    121612
    123161
    294282
    ...
    
we can also have only nonzero results (great for large mers), which prints the index, then the value

    calvin@barnabas:~/src/dna-utils$ ./kmer_total_count --nonzero -k 9 < ~/input/sample\=700013596.fa
    no input file specified with -i, reading from stdin
    0	3
    1	2
    3	3
    5	1
    ...

lastly a useful tool is having the labels generated for us, so grepping, searching and other things are easier.

    calvin@barnabas:~/src/dna/utils$ ./kmer_total_count -i ~/sample.fa -k 6 -l
    AAAAAA 552
    AAAAAC 246
    ...
    TTTTTC 102
    TTTTTG 924
    TTTTTT 4961

### kmer_counts_per_sequence



#### Licensing and Citation
Report all bugs to mutantturkey@gmail.com

Copyright 2014 Calvin Morrison, Drexel University.

If you are using any dna-utils tool for a publication
please cite your usage:

dna-utils. Drexel University, Philadelphia USA, 2014;
software available at www.github.com/EESI/dna-utils/

