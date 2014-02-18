#include <stdbool.h>

struct b_tree {
	size_t index;
	unsigned long long value;
	struct b_tree *right;
	struct b_tree *left;
};

typedef struct b_tree node;

node* search(node **tree, size_t index);
void insert(node **tree, size_t index);
void deltree(node *tree);
unsigned long long lookup(node **tree, size_t index);

void print_sparse(node *tree, bool label, bool nonzero, unsigned int kmer);
