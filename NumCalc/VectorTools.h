//
// File: VectorTools.h
// Created by: Julien Dutheil <julien.dutheil@ens-lyon.fr>
// Created on: Fri Mar 14 14:16:32 2003
//

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef _VECTORTOOLS_H_
#define _VECTORTOOLS_H_

#include "VectorExceptions.h"

#include "NumTools.h"
using namespace NumTools;

#include <vector>
using namespace std;

typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<VVdouble> VVVdouble;
typedef vector<VVVdouble> VVVVdouble;

typedef vector<int> Vint;
typedef vector<Vint> VVint;
typedef vector<VVint> VVVint;
typedef vector<VVVint> VVVVint;

namespace VectorOperators {

template<class T> vector<T> operator+ (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator+", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(int i = 0; i < size; i++) {
		result[i] = v1[i] + v2[i];
	}
	return result;
}

template<class T> vector<T> operator- (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator-", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(int i = 0; i < size; i++) {
		result[i] = v1[i] - v2[i];
	}
	return result;
}

template<class T> vector<T> operator* (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator*", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(int i = 0; i < size; i++) {
		result[i] = v1[i] * v2[i];
	}
	return result;
}

template<class T> vector<T> operator/ (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator/", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(int i = 0; i < size; i++) {
		result[i] = v1[i] / v2[i];
	}
	return result;
}



template<class T, class C> vector<T> operator+ (const vector<T> & v1, const C & c)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] + c;
	}
	return result;
}

template<class T, class C> vector<T> operator- (const vector<T> & v1, const C & c)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] - c;
	}
	return result;
}

template<class T, class C> vector<T> operator* (const vector<T> & v1, const C & c)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] * c;
	}
	return result;
}

template<class T, class C> vector<T> operator/ (const vector<T> & v1, const C & c)
{
	Vdouble result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] / c;
	}
	return result;
}



template<class T> void operator+= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] += c;
	}
}

template<class T> void operator-= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] -= c;
	}
}

template<class T> void operator*= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] *= c;
	}
}

template<class T> void operator/= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] /= c;
	}
}



template<class T, class C> void operator+= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] += c;
	}
}

template<class T, class C> void operator-= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] -= c;
	}
}

template<class T, class C> void operator*= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] *= c;
	}
}

template<class T, class C> void operator/= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] /= c;
	}
}


Vdouble operator+ (const Vdouble & v1, const Vdouble & v2) throw (DimensionException);
Vdouble operator- (const Vdouble & v1, const Vdouble & v2) throw (DimensionException);
Vdouble operator* (const Vdouble & v1, const Vdouble & v2) throw (DimensionException);
Vdouble operator/ (const Vdouble & v1, const Vdouble & v2) throw (DimensionException);

Vdouble operator+ (const Vdouble & v1, const double & c);
Vdouble operator- (const Vdouble & v1, const double & c);
Vdouble operator* (const Vdouble & v1, const double & c);
Vdouble operator/ (const Vdouble & v1, const double & c);

void operator+= (Vdouble & v1, const Vdouble & v2) throw (DimensionException);
void operator-= (Vdouble & v1, const Vdouble & v2) throw (DimensionException);
void operator*= (Vdouble & v1, const Vdouble & v2) throw (DimensionException);
void operator/= (Vdouble & v1, const Vdouble & v2) throw (DimensionException);

void operator+= (Vdouble & v1, const double & c);
void operator-= (Vdouble & v1, const double & c);
void operator*= (Vdouble & v1, const double & c);
void operator/= (Vdouble & v1, const double & c);

}

/******************************************************************************/

namespace VectorFunctions {

template<class T> T prod(const vector<T> & v1)
{
	T p = 1;
	for(unsigned int i = 0; i < v1.size(); i++) p *= v1[i];
	return p;
}
double prod(const Vdouble & v1);

template<class T> T sum(const vector<T> & v1)
{
	T p = 0;
	for(unsigned int i = 0; i < v1.size(); i++) p += v1[i];
	return p;
}
double sum(const Vdouble & v1);

template<class T> vector<double> log(const vector<T> & v1)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::log(v1[i]);
	return v2;
}
template<class T> vector<double> log(const vector<T> & v1, double base)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::log(v1[i]) / std::log(base);
	return v2;
}
Vdouble log(const Vdouble & v1);

template<class T> vector<double> log10(const vector<T> & v1)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v1.size(); i++) v2[i] = std::log10(v1[i]);
	return v2;
}
Vdouble log10(const Vdouble & v1);

template<class T> vector<T> fact(const vector<T> & v1)
{
	vector<T> v2(v1.size());
	for(unsigned int i = 0; i < v1.size(); i++) v2[i] = NumTools::fact<T>(v1[i]);
	return v2;
}

template<class T> void display(const vector<T> & v1)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		cout << v1[i] << " ";
	}
	cout << endl;
}
void display(const Vdouble & v1);

//Scalar product:
double scalar(const Vdouble & v1, const Vdouble & v2) throw (DimensionException);

//Norm of vector:
double norm(const Vdouble & v1);

//Cosinus of angle:
double cos(const Vdouble & v1, const Vdouble & v2) throw (DimensionException);

/**
 * @name Extrema.
 *
 * @{
 */
 
/**
 * @brief Template function to get the minimum value of a vector.
 *
 * The < operator must be defined for the specified class.
 *
 * @param v The input vector.
 * @return The minimum value in the vector.
 * @throw EmptyVectorException If the input vector is empty.
 */
template<class T> T min(const vector<T> & v) throw (EmptyVectorException<T>) {
	if(v.size() == 0) throw EmptyVectorException<T>("VectorFunctions::min()", & v);
	T mini = v[0];
	for(unsigned int i = 1; i < v.size(); i++)
		if(v[i] < mini) mini = v[i];
	return mini;
}

/**
 * @brief Template function to get the maximum value of a vector.
 *
 * The > operator must be defined for the specified class.
 *
 * @param v The input vector.
 * @return The maximum value in the vector.
 * @throw EmptyVectorException If the input vector is empty.
 */
template<class T> T max(const vector<T> & v) throw (EmptyVectorException<T>) {
	if(v.size() == 0) throw EmptyVectorException<T>("VectorFuntions::max()", & v);
	T maxi = v[0];
	for(unsigned int i = 1; i < v.size(); i++)
		if(v[i] > maxi) maxi = v[i];
	return maxi;
}

/**
 * @brief Template function to get the index of the minimum value of a vector.
 *
 * The < operator must be defined for the specified class.
 * The position sent is the first one matching the minimum value.
 *
 * @param v The input vector.
 * @return The position of the minimum value in the vector.
 * @throw EmptyVectorException If the input vector is empty.
 */
template<class T> int posmin(const vector<T> & v) throw (EmptyVectorException<T>) {
	T mini = min(v);
	for(unsigned int i = 0; i < v.size(); i++)
		if(v[i] == mini) return i;
	// This can't happen:
	return -1;
}

/**
 * @brief Template function to get the index of the maximum value of a vector.
 *
 * The > operator must be defined for the specified class.
 * The position sent is the first one matching the maximum value.
 *
 * @param v The input vector.
 * @return The position of the maximum value in the vector.
 * @throw EmptyVectorException If the input vector is empty.
 */
template<class T> int posmax(const vector<T> & v) throw (EmptyVectorException<T>) {
	T maxi = max(v);
	for(unsigned int i = 0; i < v.size(); i++)
		if(v[i] == maxi) return i;
	// This can't happen:
	return -1;
}

/** @} */
}

/******************************************************************************/

namespace VectorStatTools {

template<class T> double mean(const vector<T> & v1) { return sum(v1) / v1.size(); }
double mean(const Vdouble & v1);

template<class T> double center(const vector<T> & v1) { return v1 - mean(v1); }
Vdouble center(const Vdouble & v1);

double cov(const Vdouble & v1, const Vdouble & v2) throw (DimensionException);

double var(const Vdouble & v1);

double sd(const Vdouble & v1);

double cor(const Vdouble & v1, const Vdouble & v2) throw (DimensionException);

double shannon(const Vdouble & v1, double base = 2.7182818);

}

/******************************************************************************/

#endif	//_VECTORTOOLS_H_
