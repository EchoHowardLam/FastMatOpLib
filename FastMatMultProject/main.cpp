#include <iostream>

#include <algorithm>
#include <stdlib.h>
#include <time.h>

#include "complex.h"
#include "matrix.h"
#include "strassen.h"
#include "blockwise_inversion.h"
#include "polynomial_multiplication.h"

using namespace std;

#define TEST_NUM 1

void strassenTestCase(int size)
{
	long double temp[2][TEST_NUM];
	int verified = 1;
	for (int i = 0; i < TEST_NUM; i++)
	{
		long double time;
		BinaryMatrix a = Strassen::genRandomBinaryMatrix(size);
		BinaryMatrix b = Strassen::genRandomBinaryMatrix(size);

		//a.print();
		//b.print();

		time = clock();
		BinaryMatrix n = BinaryMatrix::multiplication(a, b);
		temp[0][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);
		//printf("Time for normal: %.2lf\n", (clock() - time) / ((double)CLOCKS_PER_SEC) * 1000.0);
		
		time = clock();
		BinaryMatrix s = Strassen::multiplication(a, b);
		temp[1][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);
		//printf("Time for strassen (stop at %dx%d): %.2lf\n--------------------------\n", Strassen::strassenStop, Strassen::strassenStop, (clock() - time) / ((double)CLOCKS_PER_SEC) * 1000.0);
		
		//n.print();
		//s.print();

		verified *= (BinaryMatrix::compareMatrix(n, s) == 1);
	}

	printf("Size of 2 Matrices = %d, Strassen algorithm stopping at size = %d\n", size, Strassen::strassenStop);
	printf("Normal multiplication,\n");
	for (int i = 0; i < TEST_NUM; i++)
		printf("Time: %.5lfs\n", temp[0][i]);

	printf("Strassen Algorithm,\n");
	for (int i = 0; i < TEST_NUM; i++)
		printf("Time: %.5lfs\n", temp[1][i]);

	printf("Result(s) of calculation %s\n\n", (verified == 1) ? "agree(s)" : "disagree(s)!!!");
}

void strassenTest()
{
	Strassen::strassenStop = 8;
	strassenTestCase(128);
	strassenTestCase(256);
	strassenTestCase(512);
	strassenTestCase(1024);
	Strassen::strassenStop = 16;
	strassenTestCase(128);
	strassenTestCase(256);
	strassenTestCase(512);
	strassenTestCase(1024);
	Strassen::strassenStop = 32;
	strassenTestCase(128);
	strassenTestCase(256);
	strassenTestCase(512);
	strassenTestCase(1024);
	Strassen::strassenStop = 64;
	strassenTestCase(128);
	strassenTestCase(256);
	strassenTestCase(512);
	strassenTestCase(1024);
	Strassen::strassenStop = 128;
	strassenTestCase(256);
	strassenTestCase(512);
	strassenTestCase(1024);
}

void blockInvertTestCase(int size)
{
	long double temp[2][TEST_NUM];
	int verified = 1;
	int invertible = 1;
	int success;
	for (int i = 0; i < TEST_NUM; i++)
	{
		long double time;
		BinaryMatrix a = Strassen::genRandomBinaryMatrix(size);

		time = clock();
		success = 1;
		BinaryMatrix b = BlockInvert::blockwiseInversion(a, success);
		//b.print();
		//BinaryMatrix::scalarMultiplication(b, a.det()).print();
		temp[0][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);

		time = clock();
		BinaryMatrix c = BlockInvert::blockwiseInversionStrassenMulVersion(a, success);
		temp[1][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);

		if (!success)
		{
			invertible = 0;
			continue;
		}
		verified *= (BinaryMatrix::isIdentityMatrix(Strassen::multiplication(a, b)) == 1);
		verified *= (BinaryMatrix::isIdentityMatrix(Strassen::multiplication(a, c)) == 1);
	}

	printf("Size of Matrix = %d, Blockwise inversion stopping at size = 2\n", size);
	printf("Blockwise inversion(normal multiplication),\n");
	for (int i = 0; i < TEST_NUM; i++)
		printf("Time: %.5lfs\n", temp[0][i]);

	printf("Blockwise inversion(strassen algorithm),\n");
	for (int i = 0; i < TEST_NUM; i++)
		printf("Time: %.5lfs\n", temp[1][i]);

	if (invertible)
		printf("Result(s) of calculation %s\n\n", (verified == 1) ? "agree(s)" : "disagree(s)!!!");
	else
		printf("The matrix is not invertible\n");
}

void blockInvertTest()
{
	blockInvertTestCase(2);
	blockInvertTestCase(4);
	blockInvertTestCase(8);
	blockInvertTestCase(16);
	blockInvertTestCase(32);
	blockInvertTestCase(64);
	blockInvertTestCase(128);
	blockInvertTestCase(256);
	blockInvertTestCase(512);
	blockInvertTestCase(1024);
}

void polyMulTestCase(int size)
{
	long double temp[2][TEST_NUM];
	bool verified = true;
	for (int i = 0; i < TEST_NUM; i++)
	{
		long double time;
		ValType *p = FPM::genRandomPoly(size);
		ValType *q = FPM::genRandomPoly(size);
		int newSizeNaive;
		int newSizeFast;

		//FPM::printPoly(p, size); printf("\n");
		//FPM::printPoly(q, size); printf("\n");

		time = clock();
		ComplexNumber *naiveResult = NaivePM::PolyMultiplication(p, q, size, newSizeNaive);
		temp[0][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);

		time = clock();
		ComplexNumber *fastResult = FPM::PolyMultiplication(p, q, size, newSizeFast);
		temp[1][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);

		//FPM::printPoly(naiveResult, newSizeNaive); printf("\n");
		//FPM::printPoly(fastResult, newSizeFast); printf("\n");

		// always check for the smaller ones to avoid out_of_range error, but actually newSizeFast > newSizeNaive is always true
		verified &= FPM::verifyPoly(naiveResult, fastResult, min(newSizeNaive, newSizeFast));

		delete[]p;
		delete[]q;
		delete[]naiveResult;
		delete[]fastResult;
	}

	printf("Number of Coefficient = %d\n", size);
	printf("Naive Polynomial Multiplication(n^2),\n");
	for (int i = 0; i < TEST_NUM; i++)
		printf("Time: %.5lfs\n", temp[0][i]);

	printf("Fast Polynomial Multiplication(n log n),\n");
	for (int i = 0; i < TEST_NUM; i++)
		printf("Time: %.5lfs\n", temp[1][i]);

	printf("Result(s) of calculation %s\n\n", (verified == 1) ? "agree(s)" : "disagree(s)!!!");
}

void polyMulTest()
{
	polyMulTestCase(64);
	polyMulTestCase(128);
	polyMulTestCase(256);
	polyMulTestCase(512);
	polyMulTestCase(1024);
	polyMulTestCase(2048);
	polyMulTestCase(4096);
	polyMulTestCase(8192);
	polyMulTestCase(16384);
}

int main()
{
	srand((unsigned int)time(NULL));
	//strassenTestCase(16);
	//strassenTest();
	//blockInvertTestCase(32);
	//blockInvertTest();
	//polyMulTestCase(5);
	polyMulTest();
	printf("\nAll calculations are done...\n");
	getchar();
	return 0;
}
