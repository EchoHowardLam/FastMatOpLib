#pragma once

#include <stdio.h>
#include <math.h>

#include "shared_definition.h"

class ComplexNumber
{
public:
	ValType real;
	ValType imag;
	ComplexNumber();
	ComplexNumber(ValType realVal, ValType imagVal);
	//ComplexNumber(const ComplexNumber &obj);
	ComplexNumber &operator+=(const ComplexNumber &a);
	ComplexNumber &operator-=(const ComplexNumber &a);
	ComplexNumber &operator*=(const ComplexNumber &a);
	ComplexNumber &operator*=(const ValType &a);
	ComplexNumber &operator/=(const ValType &a);
	ComplexNumber operator-() const;
	void print();
private: // unused
	ComplexNumber &operator+=(const ValType &a);
	ComplexNumber &operator-=(const ValType &a);
};

ComplexNumber operator+(ComplexNumber a, const ComplexNumber &b);
ComplexNumber operator-(ComplexNumber a, const ComplexNumber &b);

ComplexNumber operator*(ComplexNumber a, const ComplexNumber &b);
ComplexNumber &operator*=(const ValType &a, ComplexNumber b);
ComplexNumber operator*(ComplexNumber a, const ValType &b);
ComplexNumber operator*(const ValType &a, ComplexNumber b);

bool operator==(const ComplexNumber &a, const ComplexNumber &b);
bool operator!=(const ComplexNumber &a, const ComplexNumber &b);
