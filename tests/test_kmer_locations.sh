mkdir res
../kmer_locations -i test_fasta/kmer_locations.fa -k 5 -l -c  > res/kmer_locations_compliment
../kmer_locations -i test_fasta/kmer_locations.fa -k 5 -l > res/kmer_locations
../kmer_locations -i test_fasta/kmer_locations.fa -k 5 -l -c -m test_fasta/5-mer.labels.txt > res/kmer_locations_labels
