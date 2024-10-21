#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "khash.h" // khash.h from samtools/bwa/seqtk/klib
KHASH_SET_INIT_STR(s)

#define BUF_SIZE 4096

int main(int argc, char *argv[])
{
    char buf[BUF_SIZE];
    FILE *fp;
    int lineno = 0, flag = 0, ret;
    khash_t(s) *h;
    if (argc == 1) {
        fprintf(stderr, "Usage: cat in.fq | fqextract <name.lst>\n");
        return 1;
    }
    h = kh_init(s);
    fp = fopen(argv[1], "rb"); // FIXME: check fp
    while (fgets(buf, BUF_SIZE, fp))
        kh_put(s, h, strdup(buf), &ret); // FIXME: check ret
    fclose(fp);
    while (fgets(buf, BUF_SIZE, stdin)) {
        if (++lineno%4 == 1)
            flag = (kh_get(s, h, buf + 1) != kh_end(h));
        if (flag) fputs(buf, stdout);
    }
    kh_destroy(s, h); // FIXME: free keys before destroy
    return 0;
}
