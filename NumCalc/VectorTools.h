//
// File: VectorTools.h
// Created by: Julien Dutheil
// Created on: Fri Mar 14 14:16:32 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _VECTORTOOLS_H_
#define _VECTORTOOLS_H_

#include "VectorExceptions.h"

#include "NumTools.h"
using namespace NumTools;

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
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
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] / c;
	}
	return result;
}



template<class T> void operator+= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] += v2[i];
	}
}

template<class T> void operator-= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] -= v2[i];
	}
}

template<class T> void operator*= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] *= v2[i];
	}
}

template<class T> void operator/= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] /= v2[i];
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

/**
 * @brief Send the position of the first occurence of 'which'.
 *
 * Comparisons are performed using the == operator.
 * Maximum complexity: O(v.size()).
 *
 * @param v The vector to search.
 * @param which The element to search.
 * @retrun The position of which in v.
 */
template<class T> unsigned int pos(const vector<T> & v, const T & which) throw (ElementNotFoundException<T>)
{
	for(unsigned int i = 0; i < v.size(); i++)
		if(v[i] == which) return i;
	throw ElementNotFoundException<T>("VectorTools::pos.", &v, &which);
}

/**
 * @brief Send a new vector with unique elements.
 *
 * The input vector is copied, and the copy is sorted using QuickSort algorithm.
 * A one-pass loop then look for duplicates and copy unique element to a result vector.
 * The output vector is hence sorted.
 * 
 * If v is empty, it is passed 'as is' in return (after being copied).
 *
 * @param v the vector to parse.
 */
template<class T> vector<T> unique(const vector<T> & v)
{
	if(v.size() == 0) return v;
	vector<T> sortedV(v.begin(), v.end());
	sort(sortedV.begin(), sortedV.end());
	vector<T> uniq;
	uniq.push_back(sortedV[0]);
	for(unsigned int i = 1; i < sortedV.size(); i++) {
		if(sortedV[i] != sortedV[i-1]) uniq.push_back(sortedV[i]);
	}
	return uniq;
}

/**
 * @brief Tell if the vector as unique elements.
 *
 * The input vector is copied, and the copy is sorted using QuickSort algorithm.
 * A one-pass loop then look for duplicates.
 * 
 * If v is empty, the method returns 'true'.
 *
 * @param v the vector to parse.
 */
template<class T> bool isUnique(const vector<T> & v)
{
	if(v.size() == 0) return true;
	vector<T> sortedV(v.begin(), v.end());
	sort(sortedV.begin(), sortedV.end());
	for(unsigned int i = 1; i < sortedV.size(); i++) {
		if(sortedV[i] == sortedV[i-1]) return false;
	}
	return true;
}

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

template<class T> vector<double> exp(const vector<T> & v1)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::exp(v1[i]);
	return v2;
}
Vdouble exp(const Vdouble & v1);

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
