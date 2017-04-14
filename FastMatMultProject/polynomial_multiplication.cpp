#include "polynomial_multiplication.h"

ComplexNumber *FastPolynomialMultiplication::PolyMultiplication(ValType *p, ValType *q, int degree, int &new_degree)
{
	int padded_result_deg;
	ComplexNumber *pfft = FFT(p, degree, padded_result_deg);
	ComplexNumber *qfft = FFT(q, degree, padded_result_deg);
	ComplexNumber *unprocessed_r = pointwiseMultiplication(pfft, qfft, padded_result_deg);
	ComplexNumber *r = inverse_recursiveFFT(unprocessed_r, padded_result_deg);
	for (int i = 0; i < padded_result_deg; i++)
		r[i] /= (double)padded_result_deg;
	delete[]pfft;
	delete[]qfft;
	delete[]unprocessed_r;
	new_degree = padded_result_deg;
	return r;
}

ComplexNumber *FastPolynomialMultiplication::FFT(ValType *a, int len, int &new_len)
{
	if (len <= 0)
		return nullptr;
	int calLen = 1;
	while (calLen < len) calLen <<= 1;
	calLen <<= 1;
	new_len = calLen;
	ComplexNumber *t = new ComplexNumber[calLen];
	int i;
	for (i = 0; i < len; i++)
		t[i] = ComplexNumber(a[i], 0.0);
	for (; i < calLen; i++)
		t[i] = ComplexNumber();
	ComplexNumber *result = recursiveFFT(t, calLen);
	delete[]t;
	return result;
}

void FastPolynomialMultiplication::printPoly(ComplexNumber *a, int len)
{
	if (a == nullptr)
	{
		printf("Polynomial is not calculated due to error.");
		return;
	}
	for (int i = 0; i < len; i++)
	{
		if (i > 0)
			printf(" + ");
		a[i].print();
		if (i > 0)
			if (i > 1)
				printf("x^%d", i);
			else
				printf("x");

	}
	return;
}

void FastPolynomialMultiplication::printPoly(ValType *a, int len)
{
	if (a == nullptr)
	{
		printf("Polynomial is not calculated due to error.");
		return;
	}
	for (int i = 0; i < len; i++)
	{
		if (i > 0)
			printf(" + ");
		printf("%.3f", a[i]);
		if (i > 0)
			if (i > 1)
				printf("x^%d", i);
			else
				printf("x");

	}
	return;
}

ValType *FastPolynomialMultiplication::genRandomPoly(int len)
{
	ValType *o = new ValType[len];
	for (int i = 0; i < len; i++)
		o[i] = (ValType)(rand() % 1000);
	return o;
}

bool FastPolynomialMultiplication::verifyPoly(ComplexNumber *a, ComplexNumber *b, int len)
{
	for (int i = 0; i < len; i++)
		if (a[i] != b[i])
			return 0;
	return 1;
}

ComplexNumber *FastPolynomialMultiplication::pointwiseMultiplication(ComplexNumber *a, ComplexNumber *b, int len)
{
	ComplexNumber *result = new ComplexNumber[len];
	for (int i = 0; i < len; i++)
		result[i] = a[i] * b[i];
	return result;
}

ComplexNumber *FastPolynomialMultiplication::recursiveFFT(ComplexNumber *a, int len)
{
	if (len == 1)
	{
		ComplexNumber *ret = new ComplexNumber[1];
		ret[0] = a[0];
		return ret;
	}
	ComplexNumber wn(cos(2 * M_PI / len), sin(2 * M_PI / len));
	ComplexNumber w(1.0, 0.0);
	ComplexNumber *evenA = new ComplexNumber[len / 2];
	ComplexNumber *oddA = new ComplexNumber[len / 2];
	int eveni = 0, oddi = 0;
	for (int i = 0; i < len; i++)
		if (i % 2) // odd
			oddA[oddi++] = a[i];
		else // even
			evenA[eveni++] = a[i];
	ComplexNumber *evenF = recursiveFFT(evenA, len / 2);
	ComplexNumber *oddF = recursiveFFT(oddA, len / 2);
	delete[]evenA;
	delete[]oddA;
	ComplexNumber *result = new ComplexNumber[len];
	for (int k = 0; k < len / 2; k++)
	{
		result[k] = evenF[k] + w * oddF[k];
		result[k + len / 2] = evenF[k] - w * oddF[k];
		w *= wn;
	}
	delete[]evenF;
	delete[]oddF;
	return result;
}

ComplexNumber *FastPolynomialMultiplication::inverse_recursiveFFT(ComplexNumber *a, int len)
{
	if (len == 1)
	{
		ComplexNumber *ret = new ComplexNumber[1];
		ret[0] = a[0];
		return ret;
	}
	ComplexNumber wn(cos(-2 * M_PI / len), sin(-2 * M_PI / len));
	ComplexNumber w(1.0, 0.0);
	ComplexNumber *evenA = new ComplexNumber[len / 2];
	ComplexNumber *oddA = new ComplexNumber[len / 2];
	int eveni = 0, oddi = 0;
	for (int i = 0; i < len; i++)
		if (i % 2) // odd
			oddA[oddi++] = a[i];
		else // even
			evenA[eveni++] = a[i];
	ComplexNumber *evenF = inverse_recursiveFFT(evenA, len / 2);
	ComplexNumber *oddF = inverse_recursiveFFT(oddA, len / 2);
	delete[]evenA;
	delete[]oddA;
	ComplexNumber *result = new ComplexNumber[len];
	for (int k = 0; k < len / 2; k++)
	{
		result[k] = evenF[k] + w * oddF[k];
		result[k + len / 2] = evenF[k] - w * oddF[k];
		w *= wn;
	}
	delete[]evenF;
	delete[]oddF;
	return result;
}

ComplexNumber *NaivePolynomialMultiplication::PolyMultiplication(ValType *p, ValType *q, int degree, int &new_degree)
{
	new_degree = degree * 2;
	ComplexNumber *r = new ComplexNumber[new_degree]; // class declaration have already initialized and zeroed the value
	for (int i = 0; i < degree; i++)
		for (int j = 0; j < degree; j++)
			r[i + j] += ComplexNumber(p[i] * q[j], 0.0);
	return r;
}
