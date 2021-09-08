#include "../gf.h"
#include "../rng.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

poly create_polynomial(int m, int t)
{
	poly aux, polyn;
	gf r;

	int size_random = 16;
	unsigned char entropy_input[size_random];
	for (int i=0; i<size_random; i++){
        entropy_input[i] = rand();
	}
	randombytes_init(entropy_input, NULL, 256);

    unsigned char seed[1][size_random];

    // get random factors
    gf factors[t];
    int count = 0, contains = 0;
    int max_coeff_value = (1 << m);

    while (count < t) {
        do {
            randombytes(seed[0], size_random);         
            r = (gf) seed[0][0] | (gf)seed[0][1] << 8;
        } while (r > max_coeff_value || r == 0); // r < 2^m

        contains = 0;

        for (int i = 0; i < count; ++i) {
            if (r == factors[i]) {
                contains = 1;
                break;
            }
        }
        if (!contains) {
            factors[count] = r;
            ++count;
        }
    }

    // create polynomial multiplying factors
    aux.size = 2;
    aux.coeffs = malloc(sizeof(gf) * 2);
    aux.coeffs[0] = (gf) factors[0];
    aux.coeffs[1] = (gf)1;

    polyn.size = 2;
    polyn.coeffs = malloc(sizeof(gf) * 2);
    polyn.coeffs[1] = (gf)1;

    for (int i = 1; i < t; ++i) {
        polyn.coeffs[0] = (gf) factors[i];
        mult(&aux, aux, polyn);
    }
  
    return aux;
} 


int main(int argc, char** argv)
{
    if (argc != 4) {
        printf("bad args:  enter values for quantity, m, t \n");
        return 1;
    }
    if (atoi(argv[2]) > 16) {
    	printf("m must be less than 17 \n");
        return 1;
    }
    int qntdd = atoi(argv[1]);
    int m = atoi(argv[2]);
    int t = atoi(argv[3]);
    char* cm = argv[2];
    char* ct = argv[3];

    FILE *fp;

    char filename[20] = "polys-";
    char* sqntdd = argv[1];
    strcat(filename, sqntdd);
    strcat(filename, "-");
    strcat(filename, cm);
    strcat(filename, "-");
    strcat(filename, ct);
    strcat(filename, ".txt");
    fp = fopen(filename, "w+");

    fprintf(fp, "%d\n%d\n", m, t);

    for (int i = 0; i < qntdd; ++i) {
    	poly out = create_polynomial(m, t);
    	for (int i = 0; i < out.size; ++i) {
    		fprintf(fp, "%d ", out.coeffs[i]);
    	}
      fputs("\n", fp);
    }

    fclose(fp);
    return 0;
} 
