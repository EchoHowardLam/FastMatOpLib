#include "complex.h"

ComplexNumber::ComplexNumber() : real(ADD_IDENTITY), imag(ADD_IDENTITY) {}
ComplexNumber::ComplexNumber(ValType realVal, ValType imagVal) : real(realVal), imag(imagVal) {}

ComplexNumber &ComplexNumber::operator+=(const ComplexNumber &a)
{
	this->real += a.real;
	this->imag += a.imag;
	return *this;
}

ComplexNumber &ComplexNumber::operator+=(const ValType &a)
{
	this->real += a;
	return *this;
}

ComplexNumber &ComplexNumber::operator-=(const ComplexNumber &a)
{
	*this += (-a);
	return *this;
}

ComplexNumber &ComplexNumber::operator-=(const ValType &a)
{
	this->real -= a;
	return *this;
}

ComplexNumber &ComplexNumber::operator*=(const ComplexNumber &a)
{
	ValType newReal = this->real * a.real - this->imag * a.imag;
	this->imag = this->real * a.imag + this->imag * a.real;
	this->real = newReal;
	return *this;
}

ComplexNumber &ComplexNumber::operator*=(const ValType &a)
{
	this->real *= a;
	this->imag *= a;
	return *this;
}

ComplexNumber &ComplexNumber::operator/=(const ValType &a)
{
	this->real /= a;
	this->imag /= a;
	return *this;
}

ComplexNumber ComplexNumber::operator-() const
{
	return ComplexNumber(-this->real, -this->imag);
}

void ComplexNumber::print()
{
	printf("(");
	if (fabs(this->real - ADD_IDENTITY) > VERIFY_CUTOFF)
		printf("%.3f", this->real);
	if (fabs(this->imag - ADD_IDENTITY) > VERIFY_CUTOFF)
		printf("%+.3fi", this->imag);
	if (fabs(this->real - ADD_IDENTITY) <= VERIFY_CUTOFF && fabs(this->imag - ADD_IDENTITY) <= VERIFY_CUTOFF)
		printf("0");
	printf(")");
	return;
}

ComplexNumber operator+(ComplexNumber a, const ComplexNumber &b)
{
	a += b;
	return a;
}

ComplexNumber operator-(ComplexNumber a, const ComplexNumber &b)
{
	a -= b;
	return a;
}

ComplexNumber operator*(ComplexNumber a, const ComplexNumber &b)
{
	a *= b;
	return a;
}

ComplexNumber &operator*=(const ValType &a, ComplexNumber b)
{
	b *= a;
	return b;
}

ComplexNumber operator*(ComplexNumber a, const ValType &b)
{
	a *= b;
	return a;
}

ComplexNumber operator*(const ValType &a, ComplexNumber b)
{
	b *= a;
	return b;
}

bool operator==(const ComplexNumber &a, const ComplexNumber &b)
{
	if (fabs(a.real - b.real) > VERIFY_CUTOFF) return false;
	return (fabs(a.imag - b.imag) <= VERIFY_CUTOFF);
}

bool operator!=(const ComplexNumber &a, const ComplexNumber &b)
{
	if (fabs(a.real - b.real) > VERIFY_CUTOFF) return true;
	return (fabs(a.imag - b.imag) > VERIFY_CUTOFF);
}
