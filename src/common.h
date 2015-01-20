#ifndef COMMON_H
#define COMMON_H

#include "utilities.h"

/*the hash function from http://www.azillionmonkeys.com/qed/hash.html*/
uint32_t superfasthash (const char* data, int len);
#endif
