/*
  This file is for finding the roots of a polynomial
*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "params.h"
#include "gf-16.h"
#include "rng.h"
#include "amd64cpuinfo.h"
#include "djbsort/uint32_sort.h"

// Stack
//----------------------------
#define MAXSIZE STACKSIZE + STACKSIZE -1     
poly stack[MAXSIZE];     
int top = -1;            

int isempty() {

   if(top == -1)
      return 1;
   else
      return 0;
}
   
int isfull() {

   if(top == MAXSIZE)
      return 1;
   else
      return 0;
}

poly pop() {
   poly data;
	
   if(!isempty()) {
      data = stack[top];
      top = top - 1;   
      return data;
   } else {
      printf("Could not retrieve data, Stack is empty.\n");
   }
}

void push(poly data) {

   if(!isfull()) {
      top = top + 1;   
      stack[top] = data;
   } else {
      printf("Could not insert data, Stack is full.\n");
   }
}
//---------------------------- end stack


// x = (a*a) mod f
void sqr_mod(poly *x, poly a, poly f) 
{
	poly t;
	t.size = a.size + a.size-1;
	t.coeffs = malloc(sizeof(gf) * t.size);

    for (int i = 0; i < t.size; ++i) t.coeffs[i] = 0;

	gf result = 0;
	gf cof;

	for (int i = 0; i <= degree(a); ++i) {
		cof = a.coeffs[i];
		result = gf_mul(cof, cof);
		t.coeffs[2*i] = result;
		cof = result;
	}

    plainRem(x, &t, f);

    free(t.coeffs);
}

// h = tracemap of a mod f
void TraceMap(poly *h, poly a, poly f)
{
    poly res;
    poly tmp;

    p_copy(&res, &a);
    p_copy(&tmp, &a);

    long i;
    for (i = 0; i < GFBITS-1; i++) {
        sqr_mod(&tmp, tmp, f);
        GF_add(&res, res, tmp);
    }

    p_copy(h, &res);

    free(res.coeffs);
    free(tmp.coeffs);
}

int new_poly_degree_const(poly p) 
{
    int zeros = 0;
    int count = 1;
    int current;
    for (int i = degree(p); i >= 0; i--) {
       current = (!p.coeffs[i]) & count;
       zeros += current;
       count &= current;
    }
    return p.size - zeros - 1;
}

// multiply binary like integers
gf incrementFakeRoot(int n, int jump) 
{
  gf result = 0;
  while (jump != 0) {
    if (jump & 1) {
      // result ^= n;
      int a = n;
      while(a != 0) {
          int carry = result & a;
          result ^= a;
          a = carry << 1;
      }
    }

    n <<= 1;
    jump >>=1;
  }

  return result;

}

// Add diffDegree fake roots to f1, resulting in out
void addRoots(poly* out, poly f1, int expectedDegree, int diffDegree, int jump) 
{
    gf fakeroot = 1;
    poly polyn;
    polyn.size = 2;
    polyn.coeffs = malloc(sizeof(gf) * 2);

    p_copy(out, &f1);

    int adddegree = 1;
    int added = 0;
    for (int i = 1; i < expectedDegree; ++i) {
        adddegree -= !(i - (diffDegree+1));
        added += adddegree;
        fakeroot = 1 ^ (adddegree * added * jump);
        polyn.coeffs[1] = adddegree;
        polyn.coeffs[0] = fakeroot;
        mult(out, polyn, *out);
    }

    free(polyn.coeffs);
}


void iterative_probabilistic_rootfinding(uint32_t *factors, poly f1, int expectedDegree, int m) 
{
    poly h, aux, fx, trm;
    gf r;
    int roots = 0, fLoops = 0;

    int size_random = 16; // 32
    unsigned char seed[1][size_random];
    uint maxCoeff = (1 << m);


    h.coeffs = malloc(sizeof(gf) * 2);
    h.size = 2;
    h.coeffs[0] = (gf) 0;

    push(f1);

    fLoops = expectedDegree + (expectedDegree - 1);
    for (int i = 0; i < fLoops; ++i) {

        fx = pop();

        if (degree(fx) == 1) { 
            factors[roots] = fx.coeffs[0];
            ++roots;
            continue;
        }

        {
            do {
                do {
                    randombytes(seed[0], size_random);
                    r = (gf) seed[0][0] | (gf)seed[0][1] << 8;
                } while (r == 0 || r > maxCoeff);
                   
                // h = rx + 0   
                h.coeffs[1] = (gf) r;
                TraceMap(&trm, h, fx); 
                GCD(&aux, fx, trm);   

            } while (degree(aux) <= 0 || degree(aux) == degree(fx));
        }
      
        GF_div(&trm, &fx, aux);

        push(trm);
        push(aux);
   } // end for

   free(h.coeffs);
   free(trm.coeffs);
   free(aux.coeffs);
}


// x receives the roots of polynomial f
// coeffients of polynomial f are in field 2^m
void root_finding(gf *x, poly f, int expectedDegree, int m, gf fakeValue)
{    
    uint32_t factors[expectedDegree];

    for (int i = 0; i < expectedDegree; ++i) {
        factors[i] = 0;
    }

    poly newf;
    p_copy(&newf, &f);
    int degOriginal = new_poly_degree_const(newf);
    int diffDegree = expectedDegree - degOriginal;
    int jump = 63;

    // Add fake roots
    poly newf2;
    newf2.size = expectedDegree + 1;
    newf2.coeffs = malloc(sizeof(gf) * newf2.size);

    addRoots(&newf2, newf, expectedDegree, diffDegree, jump);
    newf2.size = expectedDegree + 1;
    free(newf.coeffs);


    // Find roots
    iterative_probabilistic_rootfinding(factors,newf2,expectedDegree,m);
    free(newf2.coeffs);


    // Remove fake roots from factors  
    uint32_t coeff;
    int i = 0, equal;
	gf fakeroot = jump ^ 1;
	int numberFakeRoots = diffDegree;
	int noFakeRoots = !diffDegree;
	int keepCoeff = 1;

	for (i = 0; i < expectedDegree; ++i) {
        coeff = factors[i];

        equal = !(coeff - fakeroot);
        keepCoeff = !(!noFakeRoots & equal);

        coeff = keepCoeff * coeff;
        factors[i] = coeff ^ (!keepCoeff * UINT16_MAX);

        fakeroot =  (keepCoeff * fakeroot) ^ (!keepCoeff * (1 ^((diffDegree - numberFakeRoots + 2) * jump)));
        numberFakeRoots -= !keepCoeff;
        noFakeRoots = !numberFakeRoots;
    }


	uint32_sort(factors, expectedDegree);

    keepCoeff = 1;
    for (i = 0; i < expectedDegree; ++i) {	
      	coeff = factors[i];
        keepCoeff -= !(i - degOriginal);
      	x[i] = coeff ^ (!keepCoeff * (UINT16_MAX ^ fakeValue));
    }

} // end root_finding
