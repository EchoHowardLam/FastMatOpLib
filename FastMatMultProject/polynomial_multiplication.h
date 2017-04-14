#pragma once

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "shared_definition.h"
#include "complex.h"

class FastPolynomialMultiplication
{
public:
	static ComplexNumber *PolyMultiplication(ValType *p, ValType *q, int degree, int &new_degree); // degree is the highest power of p and q
	static ComplexNumber *FFT(ValType *a, int len, int &new_len);
	static void printPoly(ComplexNumber *a, int len);
	static void printPoly(ValType *a, int len);
	static ValType *genRandomPoly(int len);
	static bool verifyPoly(ComplexNumber *a, ComplexNumber *b, int len);
private:
	static ComplexNumber *pointwiseMultiplication(ComplexNumber *a, ComplexNumber *b, int len);
	static ComplexNumber *recursiveFFT(ComplexNumber *a, int len, ComplexNumber wn, bool wn_given);
	static ComplexNumber *inverse_recursiveFFT(ComplexNumber *a, int len, ComplexNumber wn, bool wn_given);
};

typedef FastPolynomialMultiplication FPM;

class NaivePolynomialMultiplication
{
public:
	static ComplexNumber *PolyMultiplication(ValType *p, ValType *q, int degree, int &new_degree);
};

typedef NaivePolynomialMultiplication NaivePM;