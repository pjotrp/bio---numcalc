//
// File: Constraints.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Dec 25 19:35:17 2003
//

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
