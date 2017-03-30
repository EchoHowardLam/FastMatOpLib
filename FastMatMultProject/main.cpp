#include <iostream>

#include <stdlib.h>
#include <time.h>

#include "matrix.h"
#include "strassen.h"
#include "blockwise_inversion.h"

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

		//Strassen::strassenMultiplication(size, a, b);

		time = clock();
		BinaryMatrix n = BinaryMatrix::multiplication(a, b); // .print();
		temp[0][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);
		//printf("Time for normal: %.2lf\n", (clock() - time) / ((double)CLOCKS_PER_SEC) * 1000.0);
		
		time = clock();
		BinaryMatrix s = Strassen::strassenMultiplication(size, a, b); // .print();
		temp[1][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);
		//printf("Time for strassen (stop at %dx%d): %.2lf\n--------------------------\n", Strassen::strassenStop, Strassen::strassenStop, (clock() - time) / ((double)CLOCKS_PER_SEC) * 1000.0);
		
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
	long double temp[1][TEST_NUM];
	int verified = 1;
	for (int i = 0; i < TEST_NUM; i++)
	{
		long double time;
		BinaryMatrix a = Strassen::genRandomBinaryMatrix(size);

		time = clock();
		BinaryMatrix b = BlockInvert::blockwiseInversion(a);
		//b.print();
		//BinaryMatrix::scalarMultiplication(b, a.det()).print();
		temp[0][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);

		verified *= (BinaryMatrix::isIdentityMatrix(BinaryMatrix::multiplication(a, b)) == 1);
	}

	printf("Size of Matrix = %d, Blockwise inversion stopping at size = 2\n", size);
	printf("Blockwise inversion,\n");
	for (int i = 0; i < TEST_NUM; i++)
		printf("Time: %.5lfs\n", temp[0][i]);

	printf("Result(s) of calculation %s\n\n", (verified == 1) ? "agree(s)" : "disagree(s)!!!");
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

int main()
{
	srand((unsigned int)time(NULL));
	//strassenTest();
	//blockInvertTestCase(32);
	blockInvertTest();
	for (;;);
	return 0;
}