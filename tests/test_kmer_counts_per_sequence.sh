mkdir res

for i in `find input/ -type f `; do 
	../kmer_counts_per_sequence -k 6 < $i > res/`basename $i`.seq
done
 
