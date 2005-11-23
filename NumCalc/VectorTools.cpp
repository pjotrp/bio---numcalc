//
// File: VectorTools.cpp
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

// From Utils:
#include <Utils/TextTools.h>

#include "VectorTools.h"
using namespace VectorOperators;
using namespace VectorFunctions;
using namespace VectorStatTools;

// From the STL:
#include <cmath>
#include <iostream>
using namespace std;

/******************************************************************************
 * These functions are part of the VectorOperators namespace.                 *
 ******************************************************************************/

Vdouble VectorOperators::operator+ (const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator+", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	Vdouble result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] + v2[i];
	}
	return result;
}

/******************************************************************************/

Vdouble VectorOperators::operator- (const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator-", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	Vdouble result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] - v2[i];
	}
	return result;
}

/******************************************************************************/

Vdouble VectorOperators::operator* (const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator*", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	Vdouble result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] * v2[i];
	}
	return result;
}

/******************************************************************************/

Vdouble VectorOperators::operator/ (const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	unsigned int size;
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator/", v1.size(), v2.size());
	} else {
		size = v1.size();
	}
	Vdouble result(size);
	for(unsigned int i = 0; i < size; i++) {
		result[i] = v1[i] / v2[i];
	}
	return result;
}

/******************************************************************************/

Vdouble VectorOperators::operator+ (const Vdouble & v1, const double & c)
{
	Vdouble result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] + c;
	}
	return result;
}

/******************************************************************************/

Vdouble VectorOperators::operator- (const Vdouble & v1, const double & c)
{
	Vdouble result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] - c;
	}
	return result;
}

/******************************************************************************/

Vdouble VectorOperators::operator* (const Vdouble & v1, const double & c)
{
	Vdouble result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] * c;
	}
	return result;
}

/******************************************************************************/

Vdouble VectorOperators::operator/ (const Vdouble & v1, const double & c)
{
	Vdouble result(v1.size());
	for(unsigned int i = 0; i < result.size(); i++) {
		result[i] = v1[i] / c;
	}
	return result;
}

/******************************************************************************/

void VectorOperators::operator+= (Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator+=", v1.size(), v2.size());
	} 
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] += v2[i];
	}
}

/******************************************************************************/

void VectorOperators::operator-= (Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator-=", v1.size(), v2.size());
	} 
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] -= v2[i];
	}
}

/******************************************************************************/

void VectorOperators::operator*= (Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator*=", v1.size(), v2.size());
	} 
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] *= v2[i];
	}
}

/******************************************************************************/

void VectorOperators::operator/= (Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorOperators::operator/=", v1.size(), v2.size());
	} 
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] /= v2[i];
	}
}

/******************************************************************************/

void VectorOperators::operator+= (Vdouble & v1, const double & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] += c;
	}
}

/******************************************************************************/

void VectorOperators::operator-= (Vdouble & v1, const double & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] -= c;
	}
}

/******************************************************************************/

void VectorOperators::operator*= (Vdouble & v1, const double & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] *= c;
	}
}

/******************************************************************************/

void VectorOperators::operator/= (Vdouble & v1, const double & c)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		v1[i] /= c;
	}
}

/******************************************************************************/


/******************************************************************************
 *  These functions are part of the VectorFuncions namespace.                 *
 ******************************************************************************/

double VectorFunctions::prod(const Vdouble & v1)
{
	double p = 1;
	for(unsigned int i = 0; i < v1.size(); i++) p *= v1[i];
	return p;
}

/******************************************************************************/

double VectorFunctions::sum(const Vdouble & v1)
{
	double p = 0;
	for(unsigned int i = 0; i < v1.size(); i++) p += v1[i];
	return p;
}

/******************************************************************************/

Vdouble VectorFunctions::log(const Vdouble & v1)
{
	Vdouble v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::log(v1[i]);
	return v2;
}

/******************************************************************************/

Vdouble VectorFunctions::exp(const Vdouble & v1)
{
	Vdouble v2(v1.size());
	for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::exp(v1[i]);
	return v2;
}

/******************************************************************************/

Vdouble VectorFunctions::log10(const Vdouble & v1)
{
	Vdouble v2(v1.size());
	for(unsigned int i = 0; i < v1.size(); i++) v2[i] = std::log10(v1[i]);
	return v2;
}

/******************************************************************************/

void VectorFunctions::print(const Vdouble & v1, ostream & out)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		out << v1[i] << " ";
	}
	out << endl;
}

/******************************************************************************/

double VectorFunctions::scalar(const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	if(v1.size() != v2.size()) {
		throw DimensionException("VectorFunctions::scalar", v1.size(), v2.size());
	}
	double result = 0;	
	for(unsigned int i = 0; i < v1.size(); i++) {
		result += v1[i] * v2[i];
	}
	return result;
}

/******************************************************************************/

double VectorFunctions::norm(const Vdouble & v1) {
	double result = 0;
	for(unsigned int i = 0; i < v1.size(); i++) result += v1[i] * v1[i];
	return sqrt(result);
}

/******************************************************************************/

double VectorFunctions::cos(const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	return scalar(v1, v2)/(norm(v1) * norm(v2));
}

/******************************************************************************/

/******************************************************************************
 * These functions are part of the VectorStatTools namespace.                 *
 ******************************************************************************/

double VectorStatTools::mean(const Vdouble & v1) { return sum(v1) / v1.size(); }

/*****************************************************************************/

Vdouble VectorStatTools::center(const Vdouble & v1) { return v1 - mean(v1); }

/*****************************************************************************/

double VectorStatTools::cov(const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	return scalar(center(v1), center(v2)) / v1.size();
}

/*****************************************************************************/

double VectorStatTools::var(const Vdouble & v1) { return cov(v1, v1); }

/*****************************************************************************/

double VectorStatTools::sd(const Vdouble & v1) { return sqrt(var(v1)); }

/*****************************************************************************/

double VectorStatTools::cor(const Vdouble & v1, const Vdouble & v2)
throw (DimensionException)
{
	return cov(v1, v2) / ( sd(v1) * sd(v2) );
}

/*****************************************************************************/

double VectorStatTools::shannon(const Vdouble & v1, double base)
{
  double s = 0;
  for(unsigned int i = 0; i < v1.size(); i++) s += v1[i] * std::log(v1[i]) / std::log(base);
	return -s;
}

/*****************************************************************************/

