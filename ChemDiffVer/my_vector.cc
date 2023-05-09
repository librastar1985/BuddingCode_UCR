#include "my_vector.h"

using namespace std;
// #include <stdio.h>
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// #include <fstream>
// #include <cmath>

#if defined(__cplusplus)
CVector::CVector()
{
	x[0]=x[1]=x[2]=0.0;
}

CVector::CVector(double a)
{
	x[0]=x[1]=x[2] = a;
}

CVector::CVector(double a,double b, double c)
{
	x[0]=a;
	x[1]=b;
	x[2]=c;
}

double Modul(const CVector& a)
{
	return sqrt(a*a);
}

void CVector::Print()
{
	std::cout<<"("<<x[0]<<","<<x[1]<<","<<x[2]<<")\n";
}

void CVector::Printf(char *s)
{
	std::cout<<s << "(" << setiosflags(ios::fixed) << setprecision(11) <<x[0]<<", "<<x[1]<<", "<<x[2]<<")\n";
}

CVector operator+(const CVector& a, const CVector& b)
{
	return CVector(a.x[0]+b.x[0],a.x[1]+b.x[1],a.x[2]+b.x[2]);
}

CVector operator-(const CVector& a, const CVector& b)
{
	return CVector(a.x[0]-b.x[0],a.x[1]-b.x[1],a.x[2]-b.x[2]);
}

double operator*(const CVector& a, const CVector& b)
{
	return (a.x[0]*b.x[0]+a.x[1]*b.x[1]+a.x[2]*b.x[2]);
}

CVector operator*(const double& a, const CVector& b)
{
	return CVector(a*b.x[0],a*b.x[1],a*b.x[2]);
}

CVector operator*(const CVector& a, const double& b)
{
	return CVector(a.x[0]*b,a.x[1]*b,a.x[2]*b);
}

CVector operator/(const double& a, const CVector& b)
{
	return CVector(a/b.x[0],a/b.x[1],a/b.x[2]);
}

CVector operator/(const CVector& a, const double& b)
{
	return CVector(a.x[0]/b,a.x[1]/b,a.x[2]/b);
}

CVector operator>(const CVector& a, const CVector& b)
{
	return (a.x[0]>b.x[0] && a.x[1]>b.x[1] && a.x[2]>b.x[2]);
}

CVector operator<(const CVector& a, const CVector& b)
{
	return (a.x[0]<b.x[0] && a.x[1]<b.x[1] && a.x[2]<b.x[2]);
}

CVector operator>=(const CVector& a, const CVector& b)
{
	return (a.x[0]>=b.x[0] && a.x[1]>=b.x[1] && a.x[2]>=b.x[2]);
}

CVector operator<=(const CVector& a, const CVector& b)
{
	return (a.x[0]<=b.x[0] && a.x[1]<=b.x[1] && a.x[2]<=b.x[2]);
}

CVector Cross(const CVector& a, const CVector& b)
{
	// return CVector(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
	return CVector(a.x[1]*b.x[2]-a.x[2]*b.x[1],a.x[2]*b.x[0]-a.x[0]*b.x[2],a.x[0]*b.x[1]-a.x[1]*b.x[0]);
}

CVector& CVector::operator+=(CVector a)
{
    x[0] = x[0] + a.x[0];
    x[1] = x[1] + a.x[1];
    x[2] = x[2] + a.x[2];
    return *this;
}

#endif /* defined(__cplusplus) */
