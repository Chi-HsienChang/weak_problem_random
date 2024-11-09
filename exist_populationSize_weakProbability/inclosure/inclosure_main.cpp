#include "inclosure.h"

int main(int argc, char *argv[]){
    int ell = atoi(argv[1]);
    int n = atoi(argv[2]);
    int GHC = atoi(argv[3]);
    int debug = atoi(argv[4]);
    int seed = atoi(argv[5]);
    int sample_size = atoi(argv[6]);

    Inclosure A(ell, n, GHC, debug, seed, sample_size);
    A.execute();
    printf("Complete");
    return 0;
}