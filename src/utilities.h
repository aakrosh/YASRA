/*! Routines for 
 *  	Error reporting
 *  	Memory allocation
 *  	File Handling
 *  	String operations
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifndef __USE_BSD
#define __USE_BSD
#endif

#include <stdlib.h>		
#include <ctype.h>
#include <stdio.h>		/*io functions*/
#include <string.h>		/*string functions*/
#include <stdarg.h>		/*va_list and other structures*/
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>
#include <getopt.h>
#include <time.h>

#include "asserts.h"
#include "runtime.h"

#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#undef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))

# define UNUSED __attribute__((unused))

/*the memory is incremented in chunks of bytes specified by this macro*/
#define CHUNK 60

/*the return value from a routine*/
typedef enum rt_status_st {ENDFILE,SUCCESS,FAILURE}rt_status;

typedef unsigned int uint;

typedef unsigned char uchar;

typedef enum boolean{FALSE = 0, TRUE} bool;

/*the program using this header*/
extern char* argv0;

/* the time when we started the execution of this program */
extern time_t t0;

/*print the name of the program*/
void print_argv0();

/*error reporting routines*/
void fatal(const char* const msg);
void fatalf(const char* const fmt, ...);

/*memory allocation routines*/
void* ckalloc(const size_t size);
void* ckallocz(const size_t size);
void *ckrealloc(void * p, const size_t size);
void ckfree(void* const p);

/* time and memory management outputs */
void timestamp(const char* const string);
void print_usage();

/*resource allocation routine*/
void allocate_resources();

/*file handling routines*/
FILE* ckopen(const char* const name, const char* const mode);

/*string handling routines*/
char* copy_string(const char* const str);
bool same_string(const char *s, const char *t);
int compare_names(const char* const s1,const char* const s2, const int len);
char* reverse_string(char* string, const int len);
void append(uchar** parray, 	/*the pointer to the array to be used*/
            uint* const pmax, 	/*how many bytes have been allocated*/
		    uint* const plen, 	/*how many bytes have been used*/
			const int ch);		/*the byte to be added*/

/*drop in replacement for the getline function*/
signed long getline(char** lineptr, size_t* max, FILE* stream);

/*run the passed command in the shell*/
void run_command(const char* const command);
#endif
