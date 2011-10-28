/*
 * Compare our checksum list on files to the master list
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"

#define FLEN 255
#define MD5LEN 33
#define MASTER_NAMESTART 2
#define TEST_NAMESTART 27

struct filehash {
    char name[FLEN];
    char md5sum[MD5LEN];
    UT_hash_handle hh; /* makes this structure hashable */
};

FILE* file_open(char* filename) {
    FILE* fobj=fopen(filename,"r");
    if (!fobj) {
        fprintf(stderr,"Could not open file: %s\n", filename);
        exit(45);
    }
    return fobj;
}

struct filehash* load_master(char* filename) {
    char name[FLEN], md5sum[MD5LEN];
    struct filehash* files=NULL;
    int namestart=2;

    FILE* fobj=file_open(filename);

    while ((fscanf(fobj, "%s %s", md5sum, name) == 2)) {
        struct filehash* thisfile=calloc(1,sizeof(struct filehash));

        strcpy(thisfile->name, name+MASTER_NAMESTART);
        strcpy(thisfile->md5sum, md5sum);
        HASH_ADD_STR(files, name, thisfile);
    }
    fclose(fobj);

    return files;
}

struct filehash* find_file(struct filehash* files, char* name) {
    // a single file reference, don't allocate
    struct filehash* afile=NULL;
    HASH_FIND_STR(files, name, afile);
    return afile;
}

void compare_md5sums(struct filehash* master_files, char* filename) {
    struct filehash* master_file=NULL;
    char name[FLEN], md5sum[MD5LEN];

    FILE* fobj=file_open(filename);

    while ((fscanf(fobj, "%s %s", md5sum, name) == 2)) {
        master_file=find_file(master_files, name+TEST_NAMESTART);
        if (master_file==NULL) {
            fprintf(stderr,"%s not found in master list\n", name);
        } else {
            if (strcmp(md5sum, master_file->md5sum) != 0) {
                fprintf(stderr,"%s has md5sum %s instead of %s\n", 
                               name, md5sum, master_file->md5sum);
            }
        }
    }
}
int main(int argc, char** argv) {
    if (argc < 3) {
        puts("usage: checksum-compare master_list test_list\n");
        exit(45);
    }

    // this hash table will represent our files/md5sums, keyed by name
    struct filehash* master_files=NULL;

    master_files = load_master(argv[1]);
    compare_md5sums(master_files, argv[2]);

    return 0;
}
