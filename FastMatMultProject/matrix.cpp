#include "matrix.h"

BinaryMatrix::BinaryMatrix(): data(nullptr), size(0)
{
	return;
}

BinaryMatrix::BinaryMatrix(int new_size, ValType **bmat_data)
{
	data = nullptr;
	setMatrix(new_size, bmat_data);
	return;
}

BinaryMatrix::BinaryMatrix(const BinaryMatrix &obj)
{
	size = obj.size;
	data = new ValType*[size];
	for (int i = 0; i < size; i++)
	{
		data[i] = new ValType[size];
		for (int j = 0; j < size; j++)
		{
			data[i][j] = obj.data[i][j];
		}
	}
}

BinaryMatrix &BinaryMatrix::operator=(const BinaryMatrix &obj)
{
	if (this != &obj)
	{
		BinaryMatrix newObj(obj);
		std::swap(size, newObj.size);
		std::swap(data, newObj.data);
	}
	return *this;
}

BinaryMatrix::BinaryMatrix(BinaryMatrix &&obj)
{
	size = obj.size;
	data = obj.data;
	obj.data = nullptr;
}

BinaryMatrix &BinaryMatrix::operator=(BinaryMatrix &&obj)
{
	if (this != &obj)
	{
		std::swap(size, obj.size);
		std::swap(data, obj.data);
	}
	return *this;
}

BinaryMatrix::~BinaryMatrix()
{
	cleanUpMemory();
}

void BinaryMatrix::setMatrix(int new_size, ValType **bmat_data)
{
	cleanUpMemory();
	size = new_size;
	data = bmat_data;
	return;
}

ValType **BinaryMatrix::zeroMatrix(int new_size)
{
	ValType **newMat = new ValType*[new_size];
	for (int i = 0; i < new_size; i++)
	{
		newMat[i] = new ValType[new_size];
		for (int j = 0; j < new_size; j++)
		{
			newMat[i][j] = ADD_IDENTITY;
		}
	}
	return newMat;
}

ValType **BinaryMatrix::randMatrix(int new_size)
{
	ValType **newMat = new ValType*[new_size];
	for (int i = 0; i < new_size; i++)
	{
		newMat[i] = new ValType[new_size];
		for (int j = 0; j < new_size; j++)
		{
			newMat[i][j] = (ValType) (rand() % 1000) / 500;
		}
	}
	return newMat;
}

BinaryMatrix *BinaryMatrix::doubleBisection() const
{
	if (size <= 1) return nullptr;
	int newDim = size / 2;
	BinaryMatrix *ret = new BinaryMatrix[4];
	ValType **temp[4];
	for (int k = 0; k < 4; k++)
	{
		temp[k] = new ValType*[newDim];
		for (int i = 0; i < newDim; i++)
			temp[k][i] = new ValType[newDim];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i < newDim && j < newDim) {
				// top left submatrix --> [0]
				temp[0][i][j] = data[i][j];
			}
			else if (i >= newDim && j < newDim) {
				// bottom left submatrix --> [2]
				temp[2][i - newDim][j] = data[i][j];
			}
			else if (i < newDim && j >= newDim) {
				// top right submatrix --> [1]
				temp[1][i][j - newDim] = data[i][j];
			}
			else if (i >= newDim && j >= newDim) {
				// bottom right submatrix --> [3]
				temp[3][i - newDim][j - newDim] = data[i][j];
			}
		}
	}
	ret[0].setMatrix(newDim, temp[0]);
	ret[1].setMatrix(newDim, temp[1]);
	ret[2].setMatrix(newDim, temp[2]);
	ret[3].setMatrix(newDim, temp[3]);
	return ret;
}

BinaryMatrix BinaryMatrix::mergeBisection(int new_size, const BinaryMatrix *c)
{
	int oldDim = new_size / 2;
	ValType **newMat = new ValType*[new_size];
	for (int i = 0; i < new_size; i++)
	{
		newMat[i] = new ValType[new_size];
		for (int j = 0; j < new_size; j++)
		{
			if (i < oldDim && j < oldDim) {
				// top left submatrix --> [0]
				newMat[i][j] = c[0].data[i][j];
			}
			else if (i >= oldDim && j < oldDim) {
				// bottom left submatrix --> [2]
				newMat[i][j] = c[2].data[i - oldDim][j];
			}
			else if (i < oldDim && j >= oldDim) {
				// top right submatrix --> [1]
				newMat[i][j] = c[1].data[i][j - oldDim];
			}
			else if (i >= oldDim && j >= oldDim) {
				// bottom right submatrix --> [3]
				newMat[i][j] = c[3].data[i - oldDim][j - oldDim];
			}
		}
	}
	return BinaryMatrix(new_size, newMat);
}

BinaryMatrix BinaryMatrix::addition(BinaryMatrix a, const BinaryMatrix &b)
{
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			a.data[i][j] = a.data[i][j] + b.data[i][j];
	return a;
}

BinaryMatrix BinaryMatrix::subtraction(BinaryMatrix a, const BinaryMatrix &b)
{
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			a.data[i][j] = a.data[i][j] - b.data[i][j];
	return a;
}

BinaryMatrix BinaryMatrix::multiplication(const BinaryMatrix &a, const BinaryMatrix &b)
{
	BinaryMatrix newMat(a.size, zeroMatrix(a.size));
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			for (int k = 0; k < a.size; k++)
				newMat.data[i][j] += a.data[i][k] * b.data[k][j];
	return newMat;
}

BinaryMatrix BinaryMatrix::scalarMultiplication(BinaryMatrix a, const double &v)
{
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			a.data[i][j] = v * a.data[i][j];
	return a;
}

BinaryMatrix BinaryMatrix::negation(BinaryMatrix a)
{
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			a.data[i][j] = -a.data[i][j];
	return a;
}

int BinaryMatrix::compareMatrix(const BinaryMatrix &a, const BinaryMatrix &b)
{
	if (a.size != b.size) return 0;
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
			if (fabs(a.data[i][j] - b.data[i][j]) > VERIFY_CUTOFF) return 0;
	return 1;
}

int BinaryMatrix::isIdentityMatrix(const BinaryMatrix &a)
{
	for (int i = 0; i < a.size; i++)
		for (int j = 0; j < a.size; j++)
		{
			if (i == j)
			{
				if (fabs(a.data[i][j] - MULT_IDENTITY) > VERIFY_CUTOFF) return 0;
			}
			else {
				if (fabs(a.data[i][j] - ADD_IDENTITY) > VERIFY_CUTOFF) return 0;
			}
		}
	return 1;
}

void BinaryMatrix::print() const
{
	// i is rows/vertical
	// j is columns
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			//printf("%8d ", data[i][j]);
			printf("%8.3f ", data[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

ValType BinaryMatrix::det() const
{
	//if (size != 2) return MULT_IDENTITY;
	return (data[0][0] * data[1][1] - data[0][1] * data[1][0]);
}

void BinaryMatrix::cleanUpMemory()
{
	if (data == nullptr) return;
	for (int i = 0; i < size; i++)
		delete[]data[i];
	delete[]data;
	data = nullptr;
	return;
}
