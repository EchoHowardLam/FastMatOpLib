#pragma once

#include "matrix.h"
#include "strassen.h"

class BlockInvert
{
public:
	static BinaryMatrix blockwiseInversion(const BinaryMatrix &a, int &success);
	static BinaryMatrix blockwiseInversionStrassenMulVersion(const BinaryMatrix &a, int &success);
private:
	static BinaryMatrix invert2x2Matrix(BinaryMatrix a, int &success);
};
