#pragma once

#include "matrix.h"

class Strassen
{
public:
	static int strassenStop;
	static BinaryMatrix strassenMultiplication(int size, const BinaryMatrix &a, const BinaryMatrix &b);
	static BinaryMatrix *prepareCommonProcedures(int size, const BinaryMatrix &a, const BinaryMatrix &b);
	static BinaryMatrix asa(BinaryMatrix a, const BinaryMatrix &b, const BinaryMatrix &c, const BinaryMatrix &d);
	static BinaryMatrix saa(BinaryMatrix a, const BinaryMatrix &b, const BinaryMatrix &c, const BinaryMatrix &d);
	static BinaryMatrix genRandomBinaryMatrix(int size);
};
