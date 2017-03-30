#pragma once

#include "matrix.h"
#include "strassen.h"

class BlockInvert
{
public:
	static BinaryMatrix blockwiseInversion(const BinaryMatrix &a);
private:
	static BinaryMatrix invert2x2Matrix(BinaryMatrix a);
};
