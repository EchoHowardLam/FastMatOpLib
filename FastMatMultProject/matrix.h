#pragma once

#define ADD_IDENTITY 0

typedef int ValType;

class BinaryMatrix
{
public:
	BinaryMatrix();
	BinaryMatrix(int new_size, ValType **bmat_data);
	BinaryMatrix(const BinaryMatrix &obj);
	BinaryMatrix &operator=(const BinaryMatrix &obj);
	BinaryMatrix(BinaryMatrix &&obj);
	BinaryMatrix &operator=(BinaryMatrix &&obj);
	~BinaryMatrix();
	void setMatrix(int new_size, ValType **bmat_data);
	static ValType **zeroMatrix(int new_size);
	static ValType **randMatrix(int new_size);
	BinaryMatrix *doubleBisection() const;
	static BinaryMatrix mergeBisection(int new_size, const BinaryMatrix *c);
	void print() const;
	static BinaryMatrix addition(BinaryMatrix a, const BinaryMatrix &b);
	static BinaryMatrix subtraction(BinaryMatrix a, const BinaryMatrix &b);
	static BinaryMatrix multiplication(const BinaryMatrix &a, const BinaryMatrix &b);
	ValType **data;
	int size;
private:
	void cleanUpMemory();
};
