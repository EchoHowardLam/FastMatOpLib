#pragma once

#define INT_USED 0

#define ADD_IDENTITY 0.0
#define MULT_IDENTITY 1.0

typedef double ValType;

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
	static BinaryMatrix addition(BinaryMatrix a, const BinaryMatrix &b);
	static BinaryMatrix subtraction(BinaryMatrix a, const BinaryMatrix &b);
	static BinaryMatrix multiplication(const BinaryMatrix &a, const BinaryMatrix &b);
	static BinaryMatrix scalarMultiplication(BinaryMatrix a, const double &v);
	static BinaryMatrix negation(BinaryMatrix a);
	static int compareMatrix(const BinaryMatrix &a, const BinaryMatrix &b); // return 1 if a = b, otherwise return 0
	static int isIdentityMatrix(const BinaryMatrix &a); // return 1 if a = I, otherwise return 0
	void print() const;
	ValType det() const;
	ValType **data;
	int size;
private:
	void cleanUpMemory();
};
