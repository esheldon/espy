/*
 * Compare our checksum list on files to the master list
 *
 * khash is simpler but only works with basic key-value types.
 *   also declares functions for your type, which could
 *   cause symbol problems
 * uthash is a bit harder to work with but you get to work with a full struct and
 *   it is pure macro expansion.
 *
 * For my purposes this is just a key-value thing, so khash 
 * seems better.  Oddly, I see khash is slower despite indications
 * from the online benchmarks.  khash is about same as the python
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "khash.h"

#define FLEN 255
#define MD5LEN 33
#define MASTER_NAMESTART 2
#define TEST_NAMESTART 27

KHASH_MAP_INIT_STR(md5, char*)

FILE* file_open(char* filename) {
    FILE* fobj=fopen(filename,"r");
    if (!fobj) {
        fprintf(stderr,"Could not open file: %s\n", filename);
        exit(45);
    }
    return fobj;
}

khash_t(md5) * load_master(char* filename) {
    char name[FLEN], md5sum[MD5LEN];
    int namestart=2;
    int ret=0;
	khiter_t k;
    khash_t(md5) *hmap = kh_init(md5);
    FILE* fobj=file_open(filename);

    while ((fscanf(fobj, "%s %s", md5sum, name) == 2)) {
        k=kh_put(md5, hmap, strdup(name+MASTER_NAMESTART), &ret);
        // we would have to delete except we will just exit
        kh_value(hmap, k) = strdup(md5sum);
    }
    fclose(fobj);

    return hmap;
}

void compare_md5sums(khash_t(md5)* master_hash, char* filename) {
    char name[FLEN], md5sum[MD5LEN];
	khiter_t k;
    FILE* fobj=file_open(filename);

    while ((fscanf(fobj, "%s %s", md5sum, name) == 2)) {

        k=kh_get(md5, master_hash, name+TEST_NAMESTART);

        if (k == kh_end(master_hash)) {
            fprintf(stderr,"%s not found in master list\n", name);
        } else {
            char* master_md5sum = kh_value(master_hash, k);
            if (strcmp(md5sum, master_md5sum) != 0) {
                fprintf(stderr,"%s has md5sum %s instead of %s\n", 
                               name, md5sum, master_md5sum);
            }
        }
    }
}
int main(int argc, char** argv) {
    khash_t(md5) *master_hash = NULL;
    if (argc < 3) {
        puts("usage: checksum-compare master_list test_list\n");
        exit(45);
    }

    master_hash = load_master(argv[1]);
    compare_md5sums(master_hash, argv[2]);

    return 0;
}
