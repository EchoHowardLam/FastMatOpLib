#include <iostream>

#include <stdlib.h>
#include <time.h>

#include "matrix.h"
#include "strassen.h"

using namespace std;

void test(int size)
{
	long double temp[2][5];
	for (int i = 0; i < 5; i++)
	{
		long double time;
		BinaryMatrix a = Strassen::genRandomBinaryMatrix(size);
		BinaryMatrix b = Strassen::genRandomBinaryMatrix(size);

		//Strassen::strassenMultiplication(size, a, b);

		time = clock();
		BinaryMatrix::multiplication(a, b); // .print();
		temp[0][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC) * 1000.0;
		//printf("Time for normal: %.2f\n", (clock() - time) / ((double)CLOCKS_PER_SEC) * 1000.0);
		
		time = clock();
		Strassen::strassenMultiplication(size, a, b); // .print();
		temp[1][i] = (clock() - time) / ((long double)CLOCKS_PER_SEC) * 1000.0;
		//printf("Time for strassen (stop at %dx%d): %.2f\n--------------------------\n", Strassen::strassenStop, Strassen::strassenStop, (clock() - time) / ((double)CLOCKS_PER_SEC) * 1000.0);
	}

	printf("%d %d normal, \n", size, Strassen::strassenStop);
	for (int i = 0; i < 5; i++)
		printf("Time: %.2f\n", temp[0][i]);

	printf("%d %d strassen, \n", size, Strassen::strassenStop);
	for (int i = 0; i < 5; i++)
		printf("Time: %.2f\n", temp[1][i]);
}

void strassenTest()
{
	/*Strassen::strassenStop = 8;
	test(128);
	test(256);
	test(512);
	test(1024);
	Strassen::strassenStop = 16;
	test(128);
	test(256);
	test(512);
	test(1024);*/
	Strassen::strassenStop = 32;
	test(128);
	test(256);
	test(512);
	test(1024);
	Strassen::strassenStop = 64;
	test(128);
	test(256);
	test(512);
	test(1024);
	Strassen::strassenStop = 128;
	test(256);
	test(512);
	test(1024);
}

int main()
{
	srand((unsigned int)time(NULL));
	//test(1024);
	/*while (true)
	{
		BinaryMatrix *a = new BinaryMatrix[100];
		for (int i = 0; i < 100; i++)
			a[i] = BinaryMatrix(16, BinaryMatrix::zeroMatrix(16));
		delete[]a;
	}*/
	strassenTest();
	for (;;);
	return 0;
}