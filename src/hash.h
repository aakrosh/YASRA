#ifndef HASH_H
#define HASH_H

#include <assert.h>
#include "util.h"
#include "slist.h"

#define HASH_MAX_SIZE 24

/*an element in the hash table*/
struct hashEl
{
	struct hashEl* next;
	char* name;
	void* val;
};

/*abstraction of the hash table*/
struct hash
{
	struct hash* next;		/*next in the list*/
	struct hashEl** table;	/*hash buckets*/
	unsigned int mask;		/*mask hashCrc with this to get it to fit table*/
	int power_two_size;		/*size of the table in power of two*/
	int size;				/*size of the table*/
	int el_count;			/*element count in the hash table*/
};

/*abstraction to iterate through a hash table*/
struct hash_cookie
{
	struct hash* hash; 		/*hash table we are going to iterate*/
	int idx;				/*index of the current element in hash*/
	struct hashEl* next_el;  /*current element in the hash*/
};

/*return  a new hash table*/
struct hash* new_hash(const int power_of_two);


/*add a new hash value, return the hashEl*/
struct hashEl* hash_add(struct hash* const hash, const char* const name, void* const val);

/*add an integer value in hash*/
struct hashEl* hash_add_int(struct hash* const hash, const char* name, int val);

/*does this name exist in the hash table?*/
void* hash_lookup(const struct hash* const hash, const char* const name);

/*is this in the hash table, if yes return the value associated with it*/
void* hash_find_val(const struct hash* const hash, const char* const name);

/*return the interger associated with the name in an int hash*/
int hash_find_int(const struct hash* const hash, const char* const name);

/*this must be in the hash table, If not squeal and die*/
void* hash_must_find_val(const struct hash* const hash, const char* const name);

/*free up hash table*/
void free_hash(struct hash** phash);

/*return an object to be used by hash_next to traverse the table*/
struct hash_cookie hash_first(struct hash* const hash);

/*return the next entry in the hash table or NULL if no more*/
struct hashEl* hash_next(struct hash_cookie* cookie);

/*remove the node from the hash table*/
void hash_remove_el(struct hash* const hash, char* name);

/*print out the hash, for debugging purposes*/
void print_hash_el(struct hash* const hash);
#endif
