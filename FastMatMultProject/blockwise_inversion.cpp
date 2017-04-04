#include "blockwise_inversion.h"

#include <stdio.h>

BinaryMatrix BlockInvert::blockwiseInversion(const BinaryMatrix &a, int &success)
{
	if (a.size == 1)
	{
		BinaryMatrix ret = BinaryMatrix(a);
		if (ret.data[0][0] == 0)
			success = 0;
		else
			ret.data[0][0] = 1 / ret.data[0][0];
		return ret;
	}
	else if (a.size == 2) return invert2x2Matrix(a, success);
	/*
	let "i" denote inverse
	a = [X Y] = [A[0] A[1]]
		[Z W]	[A[2] A[3]]
	S = W - Z * iX * Y
	ia= [(iX + iX * Y * iS * Z * iX) (-iX * Y * iS)]
		[(-iS * Z * iX)              (iS)          ]
	*/
	int minorSize = a.size / 2;
	BinaryMatrix *A = a.doubleBisection();
	BinaryMatrix *c = new BinaryMatrix[4];
	BinaryMatrix ret;
	c[0] = blockwiseInversion(A[0], success);
	if (success)
	{
		c[1] = BinaryMatrix::negation(BinaryMatrix::multiplication(c[0], A[1]));
		c[3] = blockwiseInversion(BinaryMatrix::addition(A[3], BinaryMatrix::multiplication(A[2], c[1])), success);
		if (success)
		{
			c[2] = BinaryMatrix::negation(BinaryMatrix::multiplication(c[3], BinaryMatrix::multiplication(A[2], c[0])));
			c[0] = BinaryMatrix::addition(c[0], BinaryMatrix::multiplication(c[1], c[2]));
			c[1] = BinaryMatrix::multiplication(c[1], c[3]);
			ret = BinaryMatrix::mergeBisection(a.size, c);
		}
	}
	if (!success)
		ret = BinaryMatrix(a);
	delete[]A;
	delete[]c;
	return ret;
}

BinaryMatrix BlockInvert::blockwiseInversionStrassenMulVersion(const BinaryMatrix &a, int &success)
{
	if (a.size == 1)
	{
		BinaryMatrix ret = BinaryMatrix(a);
		if (ret.data[0][0] == 0)
			success = 0;
		else
			ret.data[0][0] = 1 / ret.data[0][0];
		return ret;
	}
	else if (a.size == 2) return invert2x2Matrix(a, success);
	/*
	let "i" denote inverse
	a = [X Y] = [A[0] A[1]]
	[Z W]	[A[2] A[3]]
	S = W - Z * iX * Y
	ia= [(iX + iX * Y * iS * Z * iX) (-iX * Y * iS)]
	[(-iS * Z * iX)              (iS)          ]
	*/
	int minorSize = a.size / 2;
	BinaryMatrix *A = a.doubleBisection();
	BinaryMatrix *c = new BinaryMatrix[4];
	BinaryMatrix ret;
	c[0] = blockwiseInversionStrassenMulVersion(A[0], success);
	if (success)
	{
		c[1] = BinaryMatrix::negation(Strassen::multiplication(c[0], A[1]));
		c[3] = blockwiseInversionStrassenMulVersion(BinaryMatrix::addition(A[3], Strassen::multiplication(A[2], c[1])), success);
		if (success)
		{
			c[2] = BinaryMatrix::negation(Strassen::multiplication(c[3], Strassen::multiplication(A[2], c[0])));
			c[0] = BinaryMatrix::addition(c[0], Strassen::multiplication(c[1], c[2]));
			c[1] = Strassen::multiplication(c[1], c[3]);
			ret = BinaryMatrix::mergeBisection(a.size, c);
		}
	}
	if (!success)
		ret = BinaryMatrix(a);
	delete[]A;
	delete[]c;
	return ret;
}

BinaryMatrix BlockInvert::invert2x2Matrix(BinaryMatrix a, int &success)
{
	if (a.size != 2) return a;
	ValType det = a.det();
	if (det == 0.0)
	{
		success = 0;
		return a;
	}
	ValType tmp = a.data[0][0];
	a.data[0][0] = a.data[1][1] / det;
	a.data[0][1] = -a.data[0][1] / det;
	a.data[1][0] = -a.data[1][0] / det;
	a.data[1][1] = tmp / det;
	return BinaryMatrix(a);
}
