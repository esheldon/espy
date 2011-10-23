#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"

#define FLEN 255
#define MD5LEN 33

struct filehash {
    char name[FLEN];
    char md5sum[MD5LEN];
    UT_hash_handle hh; /* makes this structure hashable */
};

struct filehash* load_files(char* filename, int master) {
    FILE* fobj=NULL;
    char name[FLEN], md5sum[MD5LEN];
    struct filehash* files=NULL;
    int namestart=0;

    if (master) {
        namestart=2;
    } else {
        namestart=27;
    }

    fobj=fopen(filename,"r");
    if (!fobj) {
        fprintf(stderr,"Could not open file: %s\n", filename);
    }

    while ((fscanf(fobj, "%s %s", md5sum, name) == 2)) {
        struct filehash* thisfile=calloc(1,sizeof(struct filehash));

        strcpy(thisfile->name, name+namestart);
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

void compare_md5sums(struct filehash* master_files, struct filehash* test_files) {
    // run through the test files, see if contained in master list, see
    // if md5sums are the same
    // If a file is not found, say so
    // if the md5sum is different, say so
    struct filehash* afile=NULL;
    struct filehash* tmp=NULL;
    struct filehash* master_file=NULL;

    HASH_ITER(hh, test_files, afile, tmp) {
        master_file=find_file(master_files, afile->name);
        if (master_file==NULL) {
            fprintf(stderr,"%s not found in master list\n", afile->name);
        } else {
            if (strcmp(afile->md5sum, master_file->md5sum) != 0) {
                fprintf(stderr,"%s has md5sum %s instead of %s\n", 
                               afile->name, afile->md5sum, master_file->md5sum);
            }
        }
    }
}
int main(int argc, char** argv) {
    if (argc < 3) {
        puts("usage: checksum-compare master_list test_list\n");
        exit(45);
    }

    // these will represent our files/md5sums
    struct filehash* master_files=NULL;
    struct filehash* test_files=NULL;

    //fprintf(stderr,"Loading master file list\n");
    master_files = load_files(argv[1],1);
    //fprintf(stderr,"Loading test file list\n");
    test_files  = load_files(argv[2],0);

    compare_md5sums(master_files, test_files);

    return 0;
}
