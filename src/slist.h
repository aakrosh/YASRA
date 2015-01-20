#ifndef LIST_H
#define LIST_H

#include <stdlib.h>
#include "util.h"

/*abstraction for singly linked lists*/
struct slist
{
	struct slist* next;
};

#define sladd_head(plist, node) ((node)->next = *(plist), *(plist) = node)

/*reverse order of  a list*/
void slreverse(void* plist);

/*return the number of elements in the list*/
int slcount(const void* const list);

/*remove the item from the list*/
void slremove(void* plist, void* const node);

/*return and remove the first element of the list */
void* slpop(void* plist);

/*return the last element of the list*/
void* sllast(void* list);

/*return the index'th element of the list, base 0 */
void* slelement(void* list, const int index);

/*return the index of the element from the list*/
int slindex(const void* const list, const void* const element);

/*free the list and set the pointer to the list to be NULL*/
void slfreelist(void* plist);

/* sort the linked list with qsort and a temporary array*/
void slsort(void* plist, int(*compare)(const void* const elem1, const void* const elem2));
#endif
