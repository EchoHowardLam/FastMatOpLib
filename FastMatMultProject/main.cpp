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
		BinaryMatrix c = BlockInvert::blockwiseInversion(a, success);
		temp[1][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC);

		if (!success)
		{
			invertible = 0;
			continue;
		}
		verified *= (BinaryMatrix::isIdentityMatrix(BinaryMatrix::multiplication(a, b)) == 1);
		verified *= (BinaryMatrix::isIdentityMatrix(BinaryMatrix::multiplication(a, b)) == 1);
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

int main()
{
	srand((unsigned int)time(NULL));
	//strassenTestCase(16);
	//strassenTest();
	//blockInvertTestCase(32);
	blockInvertTest();
	printf("\nAll calculations are done...\n");
	for (;;);
	return 0;
}