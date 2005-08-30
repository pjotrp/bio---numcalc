//
// File: Constraints.cpp
// Created by: Julien Dutheil
// Created on: Thu Dec 25 19:35:17 2003
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

#include "Constraints.h"

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

IncludingPositiveReal::IncludingPositiveReal(double lowerBound) {
	_lower = lowerBound;
}

bool IncludingPositiveReal::isCorrect(double value) const {
	return value >= _lower;
}

double IncludingPositiveReal::getLimit(double value) const {
	return isCorrect(value) ? value : _lower;
}

string IncludingPositiveReal::getDescription() const {
	return "[ " + TextTools::toString(_lower) + ", +inf [";
}

/******************************************************************************/

ExcludingPositiveReal::ExcludingPositiveReal(double lowerBound) {
	_lower = lowerBound;
}

bool ExcludingPositiveReal::isCorrect(double value) const {
	return value > _lower;
}

double ExcludingPositiveReal::getLimit(double value) const {
	return isCorrect(value) ? value : _lower;
}

string ExcludingPositiveReal::getDescription() const {
	return "] " + TextTools::toString(_lower) + ", +inf [";
}

/******************************************************************************/

Interval::Interval(double lowerBound, double upperBound) { 
	_lower = lowerBound;
	_upper = upperBound;
}

double Interval::getLimit(double value) const {
	if(isCorrect(value)) return value;
	else if (value <= _lower) return _lower;
	else return _upper;
}

/******************************************************************************/

IncludingInterval::IncludingInterval(double lowerBound, double upperBound):
	Interval(lowerBound, upperBound) {};

bool IncludingInterval::isCorrect(double value) const {
	return value >= _lower && value <= _upper;
}

string IncludingInterval::getDescription() const {
	return "[ " + TextTools::toString(_lower) + ", "
	            + TextTools::toString(_upper) + " ]";
}

/******************************************************************************/

ExcludingInterval::ExcludingInterval(double lowerBound, double upperBound):
	Interval(lowerBound, upperBound) {};
		
bool ExcludingInterval::isCorrect(double value) const {
	return value > _lower && value < _upper;
}

string ExcludingInterval::getDescription() const {
	return "] " + TextTools::toString(_lower) + ", "
             	+ TextTools::toString(_upper) + " [";
}

/******************************************************************************/

IncludingExcludingInterval::IncludingExcludingInterval(double lowerBound, double upperBound):
	Interval(lowerBound, upperBound) {};
		
bool IncludingExcludingInterval::isCorrect(double value) const {
	return value >= _lower && value < _upper;
}

string IncludingExcludingInterval::getDescription() const {
	return "[ " + TextTools::toString(_lower) + ", " 
	            + TextTools::toString(_upper) + " [";
}

/******************************************************************************/

ExcludingIncludingInterval::ExcludingIncludingInterval(double lowerBound, double upperBound):
	Interval(lowerBound, upperBound) {};

bool ExcludingIncludingInterval::isCorrect(double value) const {
	return value > _lower && value <= _upper;
}

string ExcludingIncludingInterval::getDescription() const {
	return "] " + TextTools::toString(_lower) + ", "
            	+ TextTools::toString(_upper) + " ]";
}

/******************************************************************************/
