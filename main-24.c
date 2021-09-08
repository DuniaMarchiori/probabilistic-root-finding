/*
   Created by Dúnia Marchiori on 09/21.
*/

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "root-24.h"
#include "gf-32.h"
#include "amd64cpuinfo.h"
#include "rng.h"


int
main(int argc, char *argv[])
{
    if (argc < 2) {
        printf("bad args: enter the polynomial  \n");
        return 1;
    }
    if (atoi(argv[argc-2]) > 32) {
        printf("m must be less than 33 \n");
        return 1;
    }

    int m = atoi(argv[argc-2]);
    int t = atoi(argv[argc-1]);
    int expectedT = t;

    gf rts[ t ];
    for (int i = 0; i < t; ++i) {
        rts[i] = 0;
    }

    int size_random = 16;
    srand( (unsigned)time(NULL) );
    unsigned char entropy_input[size_random];
    for (int i=0; i<size_random; i++){
        entropy_input[i] = rand();
    }
    randombytes_init(entropy_input, NULL, 256);

    poly locator;
    locator.size = expectedT+1;
    locator.coeffs = malloc(sizeof(gf) * locator.size);

    // set locator from argv
    for (int i = 1; i < argc-2; ++i) {
        locator.coeffs[i-1] = (gf) atoi(argv[i]);
    }

    gf fakeValue = 0;

    long long cycle_ini, cycle_end;
    cycle_ini = cpucycles_amd64cpuinfo();

    root_finding(rts, locator, expectedT, m, fakeValue);

    cycle_end = cpucycles_amd64cpuinfo();
    printf("%lld \n", cycle_end - cycle_ini);

    free(locator.coeffs);

    return 0;
}
