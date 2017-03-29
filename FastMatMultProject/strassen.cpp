#include "strassen.h"

#include <iostream>

int Strassen::strassenStop = 8;

BinaryMatrix Strassen::strassenMultiplication(int size, BinaryMatrix &a, BinaryMatrix &b)
{
	if (size == strassenStop) {
		return BinaryMatrix::multiplication(a, b);
	}
	else {
		BinaryMatrix *m = prepareCommonProcedures(size, a, b);
		BinaryMatrix *c = new BinaryMatrix[4];
		c[0] = asa(m[0], m[3], m[4], m[6]);
		c[1] = BinaryMatrix::addition(m[2], m[4]);
		c[2] = BinaryMatrix::addition(m[1], m[3]);
		c[3] = saa(m[0], m[1], m[2], m[5]);
		BinaryMatrix ret(BinaryMatrix::mergeBisection(size, c));
		delete[]m;
		delete[]c;
		return ret;
	}
}

BinaryMatrix *Strassen::prepareCommonProcedures(int size, BinaryMatrix &a, BinaryMatrix &b)
{
	BinaryMatrix *m = new BinaryMatrix[7];
	BinaryMatrix *A = a.doubleBisection();
	BinaryMatrix *B = b.doubleBisection();
	m[0] = strassenMultiplication(size / 2, BinaryMatrix::addition(A[0], A[3]), BinaryMatrix::addition(B[0], B[3]));
	m[1] = strassenMultiplication(size / 2, BinaryMatrix::addition(A[2], A[3]), B[0]);
	m[2] = strassenMultiplication(size / 2, A[0], BinaryMatrix::subtraction(B[1], B[3]));
	m[3] = strassenMultiplication(size / 2, A[3], BinaryMatrix::subtraction(B[2], B[0]));
	m[4] = strassenMultiplication(size / 2, BinaryMatrix::addition(A[0], A[1]), B[3]);
	m[5] = strassenMultiplication(size / 2, BinaryMatrix::subtraction(A[2], A[0]), BinaryMatrix::addition(B[0], B[1]));
	m[6] = strassenMultiplication(size / 2, BinaryMatrix::subtraction(A[1], A[3]), BinaryMatrix::addition(B[2], B[3]));
	delete[]A;
	delete[]B;
	return m;
}

BinaryMatrix Strassen::asa(BinaryMatrix a, const BinaryMatrix &b, const BinaryMatrix &c, const BinaryMatrix &d)
{
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			a.data[i][j] = a.data[i][j] + b.data[i][j] - c.data[i][j] + d.data[i][j];
	return a;
}

BinaryMatrix Strassen::saa(BinaryMatrix a, const BinaryMatrix &b, const BinaryMatrix &c, const BinaryMatrix &d)
{
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			a.data[i][j] = a.data[i][j] - b.data[i][j] + c.data[i][j] + d.data[i][j];
	return a;
}


BinaryMatrix Strassen::genRandomBinaryMatrix(int size)
{
	return BinaryMatrix(size, BinaryMatrix::randMatrix(size));
}