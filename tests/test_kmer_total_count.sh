
mkdir res
for i in `find input/ -type f `; do 
	../kmer_total_count -k 6 < $i > res/`basename $i`.total
done
