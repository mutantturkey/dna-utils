#include "stdio.h"
#include "limits.h"

#include "kmer_utils.h"

main() {

	int i = 0;
	int j = 0;
	unsigned long out = 0;

	for(j = 1; j < 32; j++) {
		printf("unsigned long long mer_%d(const char *seq) { return ", j);
		unsigned long multiply = 1;
		for(i = j - 1; i >= 0; i--){
			printf("(seq[%d] * %d) + ", i, multiply);
			multiply = multiply << 2;
		}
		printf(" 0; }\n");
	}

	printf("int (*return_fn())(const char * )mer_ptr(int kmer) { switch(kmer) {");
	for(j = 1; j < 32; j++)
		printf("case %d: return mer_%d;", j, j);
	printf("}");


}
