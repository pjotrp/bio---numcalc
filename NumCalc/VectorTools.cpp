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
	for(int i = 0; i < size; i++) {
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
	for(int i = 0; i < size; i++) {
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
	for(int i = 0; i < size; i++) {
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
	for(int i = 0; i < size; i++) {
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

Vdouble VectorFunctions::log10(const Vdouble & v1)
{
	Vdouble v2(v1.size());
	for(unsigned int i = 0; i < v1.size(); i++) v2[i] = std::log10(v1[i]);
	return v2;
}

/******************************************************************************/

void VectorFunctions::display(const Vdouble & v1)
{
	for(unsigned int i = 0; i < v1.size(); i++) {
		cout << v1[i] << " ";
	}
	cout << endl;
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
