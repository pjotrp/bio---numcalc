//
// File: NumTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 10 12:06:55 2003
//

#include "NumTools.h"

/******************************************************************************/
//export template<class T> T NumTools::abs(T a) { return a < 0 ? -a : a; }
double NumTools::abs(double a) { return abs<double>(a); }
float  NumTools::abs(float  a) { return abs<float >(a); }
long   NumTools::abs(long   a) { return abs<long  >(a); }
int    NumTools::abs(int    a) { return abs<int   >(a); }
short  NumTools::abs(short  a) { return abs<short >(a); }

/******************************************************************************/

//export template<class T> T NumTools::sign(T a) { return a < 0 ? -1 : (a == 0 ? 0 : 1); }
double NumTools::sign(double a) { return sign<double>(a); }
float  NumTools::sign(float  a) { return sign<float >(a); }
long   NumTools::sign(long   a) { return sign<long  >(a); }
int    NumTools::sign(int    a) { return sign<int   >(a); }
short  NumTools::sign(short  a) { return sign<short >(a); }

/******************************************************************************/

//export template<class T> T NumTools::max(T a, T b) { return a > b ? a : b; }
double NumTools::max(double a, double b) { return max<double>(a, b); }
float  NumTools::max(float  a, float  b) { return max<float >(a, b); }
long   NumTools::max(long   a, long   b) { return max<long  >(a, b); }
int    NumTools::max(int    a, int    b) { return max<int   >(a, b); }
short  NumTools::max(short  a, short  b) { return max<short >(a, b); }

/******************************************************************************/

//export template<class T> T NumTools::min(T a, T b) { return a < b ? a : b; }
double NumTools::min(double a, double b) { return min<double>(a, b); }
float  NumTools::min(float  a, float  b) { return min<float >(a, b); }
long   NumTools::min(long   a, long   b) { return min<long  >(a, b); }
int    NumTools::min(int    a, int    b) { return min<int   >(a, b); }
short  NumTools::min(short  a, short  b) { return min<short >(a, b); }

/******************************************************************************/

//export template<class T> T NumTools::sign(T a, T b) { return abs<T>(a) * sign<T>(b); }
double NumTools::sign(double a, double b) { return sign<double>(a, b); }
float  NumTools::sign(float  a, float  b) { return sign<float >(a, b); }
long   NumTools::sign(long   a, long   b) { return sign<long  >(a, b); }
int    NumTools::sign(int    a, int    b) { return sign<int   >(a, b); }
short  NumTools::sign(short  a, short  b) { return sign<short >(a, b); }

/******************************************************************************/

//export template<class T> T NumTools::sqr(T a) { return a * a; }
double NumTools::sqr(double a) { return sqr<double>(a); }
float  NumTools::sqr(float  a) { return sqr<float >(a); }
long   NumTools::sqr(long   a) { return sqr<long  >(a); }
int    NumTools::sqr(int    a) { return sqr<int   >(a); }
short  NumTools::sqr(short  a) { return sqr<short >(a); }

/******************************************************************************/

//export template<class T> void NumTools::swap(T & a, T & b) {
//	T swap = a;
//	a = b;
//	b = swap;	
//}

//export template<class T> void NumTools::shift(T & a, T & b, T c) {
//	a = b; b = c;
//}

//export template<class T> void NumTools::shift(T & a, T & b, T & c, T d) {
//	a = b; b = c; c = d;
//}

/******************************************************************************/
