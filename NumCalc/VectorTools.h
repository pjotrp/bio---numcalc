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

/**
 * @name Element-wise operations.
 */
namespace VectorOperators {

template<class T>
vector<T>	operator+ (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator+", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] + v2[i];
	}
	return result;
}

template<class T>
vector<T> operator- (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator-", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] - v2[i];
	}
	return result;
}

template<class T>
vector<T> operator* (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator*", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] * v2[i];
	}
	return result;
}

template<class T>
vector<T> operator/ (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator/", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	vector<T> result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] / v2[i];
	}
	return result;
}



template<class T, class C>
vector<T> operator+ (const vector<T> & v1, const C & c)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] + c;
	}
	return result;
}
template<class T, class C>
vector<T> operator+ (const C & c, const vector<T> & v1)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = c + v1[i];
	}
	return result;
}

template<class T, class C>
vector<T> operator- (const vector<T> & v1, const C & c)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] - c;
	}
	return result;
}
template<class T, class C>
vector<T> operator- (const C & c, const vector<T> & v1)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = c - v1[i];
	}
	return result;
}

template<class T, class C>
vector<T> operator* (const vector<T> & v1, const C & c)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] * c;
	}
	return result;
}
template<class T, class C>
vector<T> operator* (const C & c, const vector<T> & v1)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = c * v1[i];
	}
	return result;
}

template<class T, class C>
vector<T> operator/ (const vector<T> & v1, const C & c)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] / c;
	}
	return result;
}
template<class T, class C>
vector<T> operator/ (const C & c, const vector<T> & v1)
{
	vector<T> result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = c / v1[i];
	}
	return result;
}



template<class T>
void operator+= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] += v2[i];
	}
}

template<class T>
void operator-= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] -= v2[i];
	}
}

template<class T>
void operator*= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] *= v2[i];
	}
}

template<class T>
void operator/= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] /= v2[i];
	}
}



template<class T, class C>
void operator+= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] += c;
	}
}

template<class T, class C>
void operator-= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] -= c;
	}
}

template<class T, class C>
void operator*= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] *= c;
	}
}

template<class T, class C> 
void operator/= (vector<T> & v1, const C & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] /= c;
	}
}

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
 * @return The position of which in v.
 */
template<class T>
unsigned int which(const vector<T> & v, const T & which) throw (ElementNotFoundException<T>)
{
	for(unsigned int i = 0; i < v.size(); i++)
		if(v[i] == which) return i;
	throw ElementNotFoundException<T>("VectorTools::which.", &v, &which);
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
template<class T>
vector<T> unique(const vector<T> & v)
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
template<class T>
bool isUnique(const vector<T> & v)
{
	if(v.size() == 0) return true;
	vector<T> sortedV(v.begin(), v.end());
	sort(sortedV.begin(), sortedV.end());
	for(unsigned int i = 1; i < sortedV.size(); i++) {
		if(sortedV[i] == sortedV[i-1]) return false;
	}
	return true;
}

/**
 * @return The product of all elements in a vector.
 * @param v1 A vector.
 */
template<class T>
T prod(const vector<T> & v1)
{
	T p = 1;
	for(unsigned int i = 0; i < v1.size(); i++) p *= v1[i];
	return p;
}

/**
 * @return The sum of all elements in a vector.
 * @param v1 A vector.
 */
template<class T>
T sum(const vector<T> & v1)
{
	T p = 0;
	for(unsigned int i = 0; i < v1.size(); i++) p += v1[i];
	return p;
}

/**
 * @name These methods apply the corresponding function to each element
 * and return the result in a new vector.
 *
 * @{
 */
template<class T>
vector<double> log(const vector<T> & v1)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::log(v1[i]);
	return v2;
}
template<class T>
vector<double> log(const vector<T> & v1, double base)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::log(v1[i]) / std::log(base);
	return v2;
}

template<class T>
vector<double> exp(const vector<T> & v1)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::exp(v1[i]);
	return v2;
}

template<class T>
vector<double> log10(const vector<T> & v1)
{
	vector<double> v2(v1.size());
	for(unsigned int i = 0; i < v1.size(); i++) v2[i] = std::log10(v1[i]);
	return v2;
}

template<class T>
vector<T> fact(const vector<T> & v1)
{
	vector<T> v2(v1.size());
	for(unsigned int i = 0; i < v1.size(); i++) v2[i] = NumTools::fact<T>(v1[i]);
	return v2;
}
/** @} */

/**
 * @brief Print a vector to a stream.
 * @param v1 A vector.
 * @param out A stream.
 */
template<class T>
void print(const vector<T> & v1, ostream & out = cout)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		out << v1[i] << " ";
	}
	out << endl;
}

/**
 * @return The scalar product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @throw DimensionException If the two vector do not have the same length.
 */
template<class T>
T scalar(const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorFunctions::scalar", v1.size(), v2.size());
	}
	T result = 0;	
	for(unsigned int i = 0; i < v1.size(); i++) {
		result += v1[i] * v2[i];
	}
	return result;
}

/**
 * @return The scalar Kronecker product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @throw DimensionException If the two vector do not have the same length.
 */
template<class T>
vector<T> kroneckerMult(const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  unsigned int n1 = v1.size();
  unsigned int n2 = v2.size();
  vector<T> v3(n1*n2);
  for(unsigned int i = 0; i < n1; i++)
  {
    T v1i = v1[i];
    for(unsigned int j = 0; j < n2; j++)
    {
      v3[i*n2+j] = v1i*v2[j];
    }
  }
  return v3;
}

/**
 * @return The norm of a vector (\f$\sqrt{\sum_i^n x_i^2}\f$).
 * @param v1 A vector.
 */
template<class InputType, class OutputType>
OutputType norm(const vector<InputType> & v1)
{
	OutputType result = 0;
	for(unsigned int i = 0; i < v1.size(); i++)
    result += v1[i] * v1[i];
	return sqrt(result);
}
template<class InputType>
InputType norm(const vector<InputType> & v1)
{
	InputType result = 0;
	for(unsigned int i = 0; i < v1.size(); i++)
    result += v1[i] * v1[i];
	return sqrt(result);
}

/**
 * @return The cosinus of the angle of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @throw DimensionException If the two vector do not have the same length.
 */
template<class InputType, class OutputType>
OutputType cos(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
{
	return scalar<InputType,OutputType>(v1, v2)
    / (norm<InputType,OutputType>(v1) * norm<InputType,OutputType>(v2));
}
template<class InputType>
InputType cos(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
{
	return scalar<InputType>(v1, v2)
    / (norm<InputType>(v1) * norm<InputType>(v2));
}

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
template<class T>
T min(const vector<T> & v) throw (EmptyVectorException<T>)
{
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
template<class T>
T max(const vector<T> & v) throw (EmptyVectorException<T>)
{
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
template<class T>
unsigned int posmin(const vector<T> & v) throw (EmptyVectorException<T>)
{
	T mini = min(v);
	for(unsigned int i = 0; i < v.size(); i++)
		if(v[i] == mini) return i;
	// This is never reached but must be here, otherwise a warning is issued:
	return 0;
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
template<class T>
unsigned int whichmax(const vector<T> & v) throw (EmptyVectorException<T>)
{
	T maxi = max(v);
	for(unsigned int i = 0; i < v.size(); i++)
		if(v[i] == maxi) return i;
	// This is never reached but must be here, otherwise a warning is issued:
	return 0;
}

/** @} */

/**
 * @return The mean value of the vector.
 * @param v1 A vector.
 */
template<class InputType, class OutputType>
OutputType mean(const vector<InputType> & v1)
{ 
  return sum<InputType,OutputType>(v1) / (OutputType)v1.size();
}
template<class InputType>
InputType mean(const vector<InputType> & v1)
{ 
  return sum<InputType>(v1) / (InputType)v1.size();
}

/**
 * @return The median value of the vector.
 * @param v1 A vector.
 */
template<class InputType>
InputType median(const vector<InputType> & v1)
{
  InputType med = 0;
  sort(v1.begin(), v1.end());
  double size = (double)v1.size() / 2.;
  unsigned int i = v1.size() / 2;
  if (i == size)
  {
    //Vector size is pair
    med = double((v1[i-1] + v1[i]) / 2);
  }
  else
  {
    //Vector size is impair
    med = v1[i-1];
  }
  return med;
}

/**
 * @brief Set the mean of a vector to be 0.
 * 
 * @return A vector with mean 0.
 * @param v1 A vector.
 */
template<class InputType, class OutputType>
vector<OutputType> center(const vector<InputType> & v1)
{ 
  OutputType m = mean<InputType,OutputType>(v1);
  vector<OutputType> v(v1.size());
  for(unsigned int i = 0; i < v1.size(); i++)
  {
    v[i] = (OutputType)v1[i] - m;
  }
  return v;
}
template<class InputType>
vector<InputType> center(const vector<InputType> & v1)
{ 
  InputType m = mean<InputType>(v1);
  vector<InputType> v(v1.size());
  for(unsigned int i = 0; i < v1.size(); i++)
  {
    v[i] = (InputType)v1[i] - m;
  }
  return v;
}

/**
 * @return The covariance of two vectors.
 * To have a population estimate you have to multiply by \f$\frac{n}{n-1}\f$.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @throw DimensionException If the two vector do not have the same length.
 */
template<class InputType, class OutputType>
OutputType cov(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
{
	return scalar<InputType,OutputType>(
      center<InputType,OutputType>(v1),
      center<InputType,OutputType>(v2)
      ) / (OutputType)v1.size();
}
template<class InputType>
InputType cov(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
{
	return scalar<InputType>(
      center<InputType>(v1),
      center<InputType>(v2)
      ) / (InputType)v1.size();
}

/**
 * @return The variance of the vector.
 * To have a population estimate you have to multiply by \f$\frac{n}{n-1}\f$.
 * @param v1 The sample vector.
 */
template<class InputType, class OutputType>
OutputType var(const vector<InputType> & v1)
{
  return cov<InputType,OutputType>(v1, v1);
}
template<class InputType>
InputType var(const vector<InputType> & v1)
{
  return cov<InputType>(v1, v1);
}

/**
 * @return The variance of the vector.
 * To have a population estimate you have to multiply by \f$\sqrt{\frac{n}{n-1}}\f$.
 * @param v1 The sample vector.
 */
template<class InputType, class OutputType>
OutputType sd(const vector<InputType> & v1)
{
  return sqrt(var<InputType,OutputType>(v1));
}
template<class InputType>
InputType sd(const vector<InputType> & v1)
{
  return sqrt(var<InputType>(v1));
}

/**
 * @return The Pearson correlation coefficient of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @throw DimensionException If the two vector do not have the same length.
 */
template<class InputType, class OutputType>
OutputType cor(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
{
	return cov<InputType,OutputType>(v1, v2)
    / ( sd<InputType,OutputType>(v1) * sd<InputType,OutputType>(v2) );
}
template<class InputType>
InputType cor(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
{
	return cov<InputType>(v1, v2)
    / ( sd<InputType>(v1) * sd<InputType>(v2) );
}

/**
 * @return The Shannon entropy indice of the vector.
 * @param v1 The vector.
 * @param base The base of the logarithm to use.
 */
template<class T>
double shannon(const vector<T> & v1, double base = 2.7182818)
{
  T s = 0;
  for(unsigned int i = 0; i < v1.size(); i++) s += v1[i] * std::log(v1[i]) / std::log(base);
	return -s;
}

/**
 * @brief Test function.
 * @return true if all tests are passed.
 */
bool test();

}

#endif	//_VECTORTOOLS_H_

