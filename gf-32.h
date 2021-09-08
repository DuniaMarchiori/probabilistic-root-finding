/*
  This file is for polynomial functions and functions for field arithmetic
*/

#ifndef GF_32_H
#define GF_32_H

#include <stdint.h>

typedef uint32_t gf;

typedef struct {
	long size;
	gf* coeffs;
} poly;

gf gf_iszero(gf);
gf gf_add(gf, gf);
gf gf_mul(gf, gf);
gf gf_frac(gf, gf);
gf gf_inv(gf);

void p_copy(poly*, poly*);
void GF_add(poly*, poly, poly);
int degree(poly);
void mult(poly*, poly, poly);

void swap(poly*, poly*);
int IsOne(gf g);
int IsZero(poly*);

void plainRem(poly*, poly*, poly);
void GF_div(poly*, poly*, poly);
void GCD(poly*, poly, poly);

#endif

