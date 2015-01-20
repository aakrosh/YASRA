#include "slist.h"

/*reverse order of  a list*/
void slreverse(void* plist)
{
	struct slist** ppt = (struct slist**)plist;
	struct slist* newlist = NULL;
	struct slist* el = NULL;
	struct slist* next = NULL;

	next = *ppt;
	while(next != NULL){
		el = next;
		next = el->next;
		el->next = newlist;
		newlist = el;
	}
	*ppt = newlist;
}

/*return the number of elements in the list*/
int slcount(const void* const list)
{
	int count = 0;
	struct slist* iter = (struct slist*) list;

	for(; iter; ++count, iter = iter->next);
	return count;
}

/*return the last element of the list*/
void* sllast(void* list)
{	
	struct slist* node = (struct slist*)list;
	for(; node->next; node = node->next);
	return node;
}

/*return the index'th element of the list, base 0 */
void* slelement(void* list, const int index)
{
	struct slist* node = (struct slist*)list;
	int count = 0;
	
	if(index > (slcount(node)-1)){
		fatal("The index for the list cannot be greater than the number of elements in it");
	}

	while(count < index){
		node = node->next;
		count++;
	}

	return node;
}

/* return the index of the element from the list, return -1 if the element 
 * does not exist*/
int slindex(const void* const list, const void* const element)
{
	const struct slist* node = (struct slist*)list;
	int index = 0;
	while(node && (node != element)){
		node = node->next;
		index++;
	}
	
	if(node == NULL){
		return -1;
	}

	return index;
}

/*remove the item from the list*/
void slremove(void* plist, void* const node)
{
	struct slist* iter = *((struct slist**) plist);
	struct slist* t = (struct slist*) node;
	struct slist* pt = (struct slist*) node;

	if(iter == t){
		*(struct slist**)plist = t->next;
		return;
	}
	for(; iter && (iter != t); pt = iter, iter = iter->next);	
	pt->next = t->next;
	return;
}

/*return and remove the first element of the list */
void* slpop(void* plist)
{
	struct slist* node = *((struct slist**) plist);
	
	*(struct slist**)plist = node->next;
	return node;
}

/*free the list and set the pointer to the list to be NULL*/
void slfreelist(void* plist)
{
	struct slist** ppt = (struct slist**)plist;
	struct slist* next = *ppt;
	struct slist* el = NULL;

	while(next != NULL){
		el = next;
		next = el->next;
		ckfree((char*)el);
	}
	*ppt = NULL;
}

/* sort the linked list with qsort and a temporary array*/
void slsort(void* plist, int(*compare)(const void* const elem1, const void* const elem2))
{
	struct slist** pl = (struct slist**)plist;
	struct slist* list = *pl;

	int count = slcount(list);
	if(count > 1){
		struct slist* el = NULL;
		struct slist** array = NULL;
		int i = 0;

		array = ckalloc(count * sizeof(*array));
		for(el = list, i=0; el != NULL; el = el->next, i++){
			array[i] = el;
		}
		qsort(array, count, sizeof(array[0]), compare);
		list = NULL;
		for( i = 0; i < count; i++){
			array[i]->next = list;
			list = array[i];
		}
		ckfree(array);
		slreverse(&list);
		*pl = list;
	}
}
