/* hash a string key. The code code is mostly borrowed from Jim kents library 
 * at Santa Cruz
 */

#include "hash.h"

/*return a  new hash table for out purposes*/
struct hash* new_hash(int power_two_size)
{
	struct hash* hash = ckallocz(sizeof(struct hash));
	if(0 == power_two_size){
		power_two_size = 12;
	}
	assert(power_two_size < HASH_MAX_SIZE && power_two_size > 0);
	hash->power_two_size  = power_two_size;
	hash->size = (1<<power_two_size);
	hash->mask = hash->size - 1;
	hash->table = ckallocz(hash->size * sizeof(struct hashEl));
	return hash;
}

/*compute the hash value of a string*/
static unsigned int hash_string(const char* const string)
{
	const char *keyStr = string;
	unsigned int result = 0;
	int c;

	while ((c = *keyStr++) != '\0'){
 	   result += (result<<3) + c;
    }
	return result;
}

/*add an element to the hash table (generic function)*/
static struct hashEl* generic_hash_add(struct hash* const hash, const char* const name, const int name_size, void* const val)
{
	struct hashEl* el = ckalloc(sizeof(struct hashEl));
	int hashVal = (hash_string(name) & hash->mask);
	el->name = ckallocz((strlen(name)+1)*sizeof(char));
	memcpy(el->name, name, name_size);
	el->val = val;
	el->next = hash->table[hashVal];
	hash->table[hashVal] = el;
	hash->el_count++;
	return el;
}

/*add an element to the hash table*/
struct hashEl* hash_add(struct hash* const hash, const char* name, void* const val)
{
	return generic_hash_add(hash, name, strlen(name), val); 
}

/*add an integer value in hash*/
struct hashEl* hash_add_int(struct hash* const hash, const char* name, int val)
{
	char* pt = NULL;
	return hash_add(hash, name, pt+val);
}

/*does this name exist in the hash table?*/
void* hash_lookup(const struct hash* const hash, const char* const name)
{
	struct hashEl* hel = hash->table[hash_string(name)&hash->mask];
	while(hel != NULL){
		if(same_string(hel->name, name)){
			break;
		}
		hel = hel->next;
	}
	return hel;
}

/*return NULL or the value associated with this string from the hash table if it
 * exits*/
void* hash_find_val(const struct hash* const hash, const char* const name)
{
	struct hashEl* hel = hash_lookup(hash, name);
	if(hel == NULL){
		return NULL;
	}

	return hel->val;
}

/*return the interger associated with the name in an int hash*/
int hash_find_int(const struct hash* const hash, const char* const name)
{
	struct hashEl *hel = hash_lookup(hash, name);
	if(hel == NULL){
	    return INT_MAX;
	}

	char* a = NULL;
	return (char*)hel->val - a;
}

/*this must be in the hash table, If not squeal and die*/
void* hash_must_find_val(const struct hash* const hash, const char* const name)
{
	struct hashEl* hel = hash_lookup(hash, name);
	if(hel == NULL){
		fatalf("Did not find %s in the hash", name);
		return NULL;
	}

	return hel->val;
}



/*return an object to be used by hash_next to traverse the table*/
struct hash_cookie hash_first(struct hash* const hash)
{
	struct hash_cookie cookie;
	cookie.hash = hash;
	cookie.idx = 0;
	cookie.next_el = NULL;

	/*find first entry*/
	for(cookie.idx = 0; cookie.idx < hash->size && hash->table[cookie.idx] == NULL; cookie.idx++){
		continue;

	}
	if(cookie.idx < hash->size){
		cookie.next_el = hash->table[cookie.idx];
	}
	return cookie;
}

/*return the next entry in the hash table or NULL if no more*/
struct hashEl* hash_next(struct hash_cookie* cookie)
{
	struct hashEl* ret = cookie->next_el;
	if(ret == NULL){
		return NULL;
	}

	/*find next entry*/
	cookie->next_el = ret->next;
	if(cookie->next_el == NULL){
		for(cookie->idx++; (cookie->idx < cookie->hash->size) && (cookie->hash->table[cookie->idx] == NULL); cookie->idx++){
			continue;
		}
		if(cookie->idx < cookie->hash->size){
			cookie->next_el = cookie->hash->table[cookie->idx];
		}
	}
	return ret;
}

/*free up the hash table*/
void free_hash(struct hash**  phash)
{	
	struct hash* hash = *phash;
	struct hash_cookie cookie = hash_first(hash);
	struct hashEl* temp = NULL;
	struct hashEl* el = hash_next(&cookie);
	while(el){
		temp = hash_next(&cookie);
		ckfree(el->name);
		ckfree(el);
		el = temp;
	}
	
	ckfree(hash->table);
	ckfree(hash);
	*phash = NULL;
}

/*print out the hash names, for debugging purposes only*/
void print_hash_el(struct hash* const hash)
{
	struct hash_cookie cookie = hash_first(hash);
	struct hashEl* el = hash_next(&cookie);
	while(el){
		printf("%s\n", el->name);
		el = hash_next(&cookie);
	}
}

/*remove the element from the hash*/
void hash_remove_el(struct hash* const hash, char* name)
{
	struct hashEl* hel = NULL;
	struct hashEl** pbucket = &hash->table[hash_string(name)&hash->mask];
	for(hel = *pbucket; hel != NULL; hel = hel->next){
		if(same_string(hel->name, name)){
			break;
		}
	}

	if(hel == NULL){
		fatalf("Trying to remove non-existant %s from hash", name);
	}
	slremove(pbucket, hel);
	ckfree(hel->name);
	ckfree(hel);
	hash->el_count--;
	return;
}
