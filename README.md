# Towards constant-time probabilistic root finding for code-based cryptography


## Overview

This repository implements a proposal for a constant-time probabilistic root finding algorithm. It is a work in progress and it is not recommended to be used in a production environment. 


## Use

After cloning the repository, build using `./build`. It is important to change the values in `param.h` when necessary before building.
The `main` executable created expects the following parameters:

`./main <polynomial> <m> <expected_degree>`
- `<polynomial>`: the input polynomial
- `<m>`: the field size of the coefficients of the input polynomial is 2^m
- `<expected_degree>`: the polynomial degree expected by the cryptosystem

The output is the number of CPU cycles that the root finding algorithm took.


### Polynomials

The input polynomials is a string of coefficients, in the form `c_0 c_1 c_2 ... c_n` for a polynomial of degree `n` . Each `c_i` is the integer value of the binary polynomial in GF(2^m).
In `polynomials/` it is possible to create files with random polynomials with the script in `create_polys.sh`. The parameters are the following:

`sh create_polys.sh <quantity> <m> <begin> <end>`
- `quantity>`: the number of polynomials per file (one file for each degree)
- `<m>`: the field size of the coefficients of the polynomials is 2^m
- `<begin>`: the value which the degree range starts
- `<end>`: the value which the degree range stops (inclusive)


### Field size

Field sizes up to 2^24 are supported, but only 2^m, m =[13,16,18,20] are implemented. In order to use bigger sizes (up to 2^32), files `gf-size.c` and `root-32.c` must be created. In `root-32.c`it is important to change the creation of the random coefficient and in `gf-size.c` it is important to adapt methods `gf_mul` and `gf_frac`.
