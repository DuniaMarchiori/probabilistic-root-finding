/*
  This file is for polynomial functions and functions for field arithmetic
  Based on Classic McEliece's NIST submission code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "gf-32.h"
#include "params.h"
#include <string.h>

#define MIN(a,b) (((a)<(b))?(a):(b))


gf gf_iszero(gf a)
{
	uint64_t t = a;

	t -= 1;
	t >>= 35;

	return (gf) t;
}

gf gf_add(gf in0, gf in1)
{
	return in0 ^ in1;
}

gf gf_mul(gf in0, gf in1)
{
	int i;

	uint64_t tmp;
	uint64_t t0;
	uint64_t t1;
	uint64_t t;

	t0 = in0;
	t1 = in1;

	tmp = t0 * (t1 & 1);

	for (i = 1; i < GFBITS; i++)
		tmp ^= (t0 * (t1 & (1 << i)));

	t = tmp & 0x7FFF000000;
	tmp ^= (t >> 17) ^ (t >> 20);

	t = tmp & 0x0000F00000;
	tmp ^= (t >> 17) ^ (t >> 20);

	return tmp & GFMASK;
}

/* input: field element in */
/* return: (in^2)^2 */
static inline gf gf_sq2(gf in)
{
	gf out = gf_mul(in,in);
	return gf_mul(out,out);
}

/* input: field element in, m */
/* return: (in^2)*m */
static inline gf gf_sqmul(gf in, gf m)
{
	gf out = gf_mul(in,in);
	return gf_mul(out,m);
}

/* input: field element in, m */
/* return: ((in^2)^2)*m */
static inline gf gf_sq2mul(gf in, gf m)
{
	gf out = gf_sq2(in);
	return gf_mul(out,m);
}

/* input: field element den, num */
/* return: (num/den) */
gf gf_frac(gf den, gf num)
{
	gf tmp_11;
	gf tmp_1111;
	gf out;

	tmp_11 = gf_sqmul(den, den); // ^11
	tmp_1111 = gf_sq2mul(tmp_11, tmp_11); // ^1111
	out = gf_sq2(tmp_1111); 
	out = gf_sq2mul(out, tmp_1111); // ^11111111
	out = gf_sq2(out);
	out = gf_sq2mul(out, tmp_1111); // ^111111111111

	out = gf_sq2(out);
	out = gf_sq2mul(out, tmp_1111);  // 16

	out = gf_sq2mul(out, tmp_11); //18

	out = gf_mul(out, out);
	out = gf_mul(out, den); //19

	return gf_sqmul(out, num); // ^1111111111110 = ^-1
}

gf gf_inv(gf den)
{
	return gf_frac(den, ((gf) 1));
}

int degree(poly g) {
	return g.size - 1;
}

void update_poly_degree_const(poly* p) {
    int zeros = 0;
    int count = 1;
    int current;
    for (int i = degree(*p); i >= 0; i--) {
       current = (!p->coeffs[i]) & count;
       zeros += current;
       count &= current;
    }
    p->size -= zeros;
}

void p_copy(poly *out, poly *p)
{
    out->size = p->size;
    if (out->size >= 0) {

    	out->coeffs = malloc(sizeof(gf) * out->size);
	    memcpy(out->coeffs, p->coeffs, (sizeof(gf) * (p->size)));
    }
}

void GF_add(poly *out, poly in0, poly in1)
{
	poly sum;

	if (IsZero(&in0)) {
		p_copy(out, &in1);
		return;
	}
	else if (IsZero(&in1)) {
		p_copy(out, &in0);
		return;
	}

	if (degree(in0) > degree(in1))
	{
		p_copy(&sum, &in0);
		for (int i = 0; i <= degree(in1); ++i)
        {
        	sum.coeffs[i] ^= in1.coeffs[i];
        }
	} else {
		p_copy(&sum, &in1);
		for (int i = 0; i <= degree(in0); ++i)
        {
        	sum.coeffs[i] ^= in0.coeffs[i];
        }
	}

	p_copy(out, &sum);

	free(sum.coeffs);
}


void mult(poly *out, poly in0, poly in1)
{
	int i = 0, j = 0;

	poly mul;
	mul.size = in0.size  + in1.size - 1;
	if (mul.size < 0) mul.size = 0;
	mul.coeffs = malloc(sizeof(gf) * mul.size);

	for (i = 0; i < mul.size; ++i) mul.coeffs[i] = 0;

	for (i=0; i <= degree(in0); i++) 
   { 
		for (j = 0; j <= degree(in1); ++j)
	    {
	    	mul.coeffs[i+j] ^= gf_mul(in0.coeffs[i], in1.coeffs[j]);
	    }
	}

	p_copy(out, &mul);
	free(mul.coeffs);
}

int IsZero(poly *g) 
{
	return g->size == 0;
}

int IsOne(gf g) {
	return (g == 1);
}

void swap(poly *a, poly *b) 
{
	poly temp = *a;
	*a = *b;
	*b = temp;
}

void truncate(poly* f, poly a, int m) 
{
	p_copy(f, &a);
	int size = a.size;
    f->size = MIN(size, m);
}

void reverse(poly* f, poly a, int d) 
{
   int i, j, n, m;

   n = d+1;
   m = a.size;

   f->size = n;
   f->coeffs = malloc(sizeof(gf) * f->size);

   for (i = 0; i < n; i++) {
      j = d-i;
      if (j < 0 || j >= m)
         f->coeffs[i] = 0;
      else
         f->coeffs[i] = a.coeffs[j];
   }
}

int divstepsx(poly *f, poly *g, int n, int t, int initialDelta, poly initialF, poly initialG) 
{
	int delta = initialDelta;
    poly aux, aux2;

    truncate(f, initialF, t);
    truncate(g, initialG, t);

	poly x;
	x.size = 2;
	x.coeffs = malloc(sizeof(gf) * 2);
	x.coeffs[0] = (gf) 0;
	x.coeffs[1] = (gf) 1;
	poly f0, g0;
	f0.size = 1;
	f0.coeffs = malloc(sizeof(gf));
	g0.size = 1;
	g0.coeffs = malloc(sizeof(gf));
   


	while (n > 0) {
		truncate(f, *f, t);

		if (delta > 0 && degree(*g) >= 0 && !gf_iszero(g->coeffs[0])) {
			delta = - delta;
			p_copy(&aux, f);
			p_copy(f, g);
			p_copy(g, &aux);
    	} // end if


	    f0.coeffs[0] = f->coeffs[0];
	    g0.coeffs[0] = g->coeffs[0];
	    delta = 1 + delta;

	    mult(&aux, f0, *g);
	    mult(&aux2, g0, *f);
	    GF_add(&aux, aux, aux2);

	    GF_div(g, &aux, x);

	    n = n-1;
	    t = t-1;
	    truncate(g, *g, t);
  	} // end while

	free(aux.coeffs);
	free(aux2.coeffs);
	free(x.coeffs);
	free(f0.coeffs);
	free(g0.coeffs);

  	return delta;
}

void GCD(poly *d, poly u, poly v)
{
	poly u1, v1;

	p_copy(&u1, &u);
	p_copy(&v1, &v);

	if (degree(u1) == degree(v1)) {
		if (IsZero(&u1)) {
			return;
		}

	  plainRem(&v1, &v1, u1);
	} else if (degree(u1) < degree(v1)) {
		swap(&u1, &v1);
	}


	int deg = degree(u1);

	poly f, g;
	reverse(&f, u1, deg);
	reverse(&g, v1, deg-1);

	int delta = divstepsx(&f, &g, 2*deg-1, 3*deg-1, 1, f, g);

	poly f0;
	f0.size = 1;
	f0.coeffs = malloc(sizeof(gf));
	f0.coeffs[0] = f.coeffs[0];

	reverse(&g, f, delta/2);
	GF_div(d, &g, f0);


	free(u1.coeffs);
	free(v1.coeffs);
	free(f0.coeffs);
	free(f.coeffs);
	free(g.coeffs);
}

void GF_div(poly *quo, poly *p, poly d)
/* p: poly;  d: divisor;  r: remainder; returns quotient quo */
{
    int i, j;
    int power = degree(*p) - degree(d);
    gf ratio;

    if (power < 0) {
    	quo->size = 1;
    	quo->coeffs = malloc(sizeof(gf));
    	quo->coeffs[0] = (gf) 0;
    	printf("sai div \n");
    	return;
    }

    poly q, r;
	q.size = power+1;
	q.coeffs = malloc(sizeof(gf) * q.size);
    p_copy(&r, p);

    gf d_inv = gf_inv(d.coeffs[degree(d)]);

    for (i = degree(*p); i >= degree(d); i--) {
    		q.coeffs[i - degree(d)] = ratio = gf_mul(r.coeffs[i], d_inv);
            r.coeffs[i] = 0;

            for (j = 0; j < degree(d); j++) {
            	r.coeffs[i - degree(d) + j] ^= gf_mul(d.coeffs[j], ratio);
            }
    }

    p_copy(quo, &q);

    free(r.coeffs);
    free(q.coeffs);
}

void plainRem(poly *remp, poly *p, poly d) 
{
	int LCIsOne;
	long da, db, i, j;
	gf LCInv, t;

	da = degree(*p);
    db = degree(d);

    if (da < db) {
		p_copy(remp, p);
		return;
    }

    LCIsOne = IsOne(d.coeffs[db]);

	if (!LCIsOne) {
    	LCInv = gf_inv(d.coeffs[db]);
   	}

	poly x;
	p_copy(&x, p);


	for (i = da - db; i >= 0; i--) {
		t = x.coeffs[i+db];
		if (!LCIsOne) {
			t = gf_mul(t, LCInv);
		}

		for (j = db-1; j >= 0; j--) {
			x.coeffs[i+j] ^= gf_mul(t, d.coeffs[j]);;
		}
	}

    remp->size = db;
    remp->coeffs = malloc(sizeof(gf) * remp->size);
    for (i = 0; i < remp->size; i++) {
    	remp->coeffs[i] = x.coeffs[i];
    }

    update_poly_degree_const(remp);

    if (remp->size < 0) {
    	remp->size = 0;
    	remp->coeffs = malloc(sizeof(gf));
    	return;
    }

    free(x.coeffs);
}
