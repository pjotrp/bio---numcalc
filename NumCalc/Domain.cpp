//
// File: Domain.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Feb  4 19:12:44 2004
//

#include "Domain.h"


Domain::Domain(double a, double b, unsigned int n)
{
	double mini = min(a, b);
	double maxi = max(a, b);
	double w = (maxi - mini) / n;
	_bounds.resize(n + 1);
	_bounds[0] = mini;
	_midPoints.resize(n);
	for(unsigned int i = 1; i < n + 1; i++) {
		_bounds[i] = mini + i * w;
		_midPoints[i - 1] = mini + (i - .5) * w;
	}
}

Domain::Domain(Vdouble bounds): _bounds(bounds)
{
    unsigned int n = _bounds.size();
    _midPoints.resize(n-1);
    for(unsigned int i = 0; i < n-1; i++) {
        _midPoints[i] = (_bounds[i] + _bounds[i+1]) / 2.;
    }
}

Domain::Domain(Vdouble bounds, Vdouble midPoints) throw (Exception): _bounds(bounds)
{
    if(bounds.size() != midPoints.size() + 1) throw Exception("Domain::Domain(). Number of midpoints must equal number of bounds - 1.");
    _midPoints = midPoints;
    // Check if midPoint are really midPoints ;-)
    for(unsigned int i = 0; i < _midPoints.size(); i++) {
        if(!(_midPoints[i] >= _bounds[i] && _midPoints[i] < _bounds[i+1]))
            throw Exception("Domain::Domain(). Midpoint " +
                TextTools::toString((int)i) +
                " = " +
                TextTools::toString(_midPoints[i]) +
                " does not belong to interval [" +
                TextTools::toString(_bounds[i]) +
                ", " +
                TextTools::toString(_bounds[i+1]) +
                "]."
            );
    }
}

double Domain::getLowerBound() const { return * _bounds.begin(); }

double Domain::getUpperBound() const { return * _bounds.rbegin(); }

double Domain::getLowerValue() const { return * _midPoints.begin(); }

double Domain::getUpperValue() const { return * _midPoints.rbegin(); }

unsigned int Domain::getSize() const { return _midPoints.size(); }

double Domain::getBound(unsigned int i) const { return _bounds[i]; }

double Domain::getValue(unsigned int i) const { return _midPoints[i]; } 

double Domain::getNearestValue(double x) const throw (OutOfRangeException) 
{
	if(x < getLowerBound() || x >= getUpperBound())
		throw OutOfRangeException("Domain::getNearestValue", x, getLowerBound(), getUpperBound());
	for(unsigned int i = 1; i < _bounds.size(); i++)
		if(x < _bounds[i])
			return _midPoints[i - 1];
	// This line can't be reached:
	return 0;
}

int Domain::getIndex(double x) const throw (OutOfRangeException)
{
	if(x < getLowerBound() || x >= getUpperBound())
		throw OutOfRangeException("Domain::getIndex", x, getLowerBound(), getUpperBound());
	for(unsigned int i = 1; i < _bounds.size(); i++)
		if(x < _bounds[i])
			return i - 1;
	// This line can't be reached:
	return -1;
}
