#if !defined(_VECTOR_H)
#define _VECTOR_H

// #include <gdecs/gdecs.h>
#include "cdecs.h"

// #include <vmalloc.h>
/**
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
**/

// #include <stdio.h>
// #include <math.h>
// #include <cstdlib>
// #include <iostream>
// #include <fstream>
// #include <cmath>
// #define MAX(_a,_b) ((_a)>(_b) ? (_a) : (_b))  

// #define MIN(_a,_b) ((_a)<(_b) ? (_a) : (_b))

#if defined(__cplusplus)
class CVector
{
public:
	// double x,y,z;
	double x[3];
	CVector();
	CVector(double a);
	CVector(double a,double b, double c);
	double GetX(){ return x[0];}
	double GetY(){ return x[1];}
	double GetZ(){ return x[2];}
	void   Printf(char *name);

	friend CVector operator+(const CVector& a, const CVector& b);
	friend CVector operator-(const CVector& a, const CVector& b);
	friend double operator*(const CVector& a, const CVector& b);
	friend CVector operator*(const double& a, const CVector& b);
	friend CVector operator*(const CVector& a, const double& b);
	friend CVector operator/(const double& a, const CVector& b);
	friend CVector operator/(const CVector& a, const double& b);
	friend CVector operator<(const CVector& a, const CVector& b);
	friend CVector operator>(const CVector& a, const CVector& b);	
	friend CVector operator<=(const CVector& a, const CVector& b);
	friend CVector operator>=(const CVector& a, const CVector& b);	
	friend CVector Cross(const CVector& a, const CVector& b);
	friend double Modul(const CVector& a);

        CVector& operator+= (CVector);

	void Print();
};
#endif /* defined(__cplusplus) */

#endif /* if !defined(_VECTOR_H) */
