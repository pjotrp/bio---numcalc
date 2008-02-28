//
// File: Domain.cpp
// Created by: Julien Dutheil
// Created on: Wed Feb  4 19:12:44 2004
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

#include "Domain.h"

// From Utils:
#include <Utils/TextTools.h>

using namespace bpp;

Domain::Domain(double a, double b, unsigned int n) throw (Exception)
{
  if(n == 0) throw Exception("Domain::constructor1. Number of classes should be > 0.");
	double mini = min(a, b);
	double maxi = max(a, b);
	double w = (maxi - mini) / n;
	_bounds.resize(n + 1);
	_bounds[0] = mini;
	_midPoints.resize(n);
	for(unsigned int i = 1; i < n + 1; i++)
  {
		_bounds[i] = mini + i * w;
		_midPoints[i - 1] = mini + (i - .5) * w;
	}
}

Domain::Domain(Vdouble bounds) throw (Exception): _bounds(bounds)
{
  unsigned int n = _bounds.size();
  _midPoints.resize(n-1);
  for(unsigned int i = 0; i < n-1; i++) {
		if(bounds[i+1] <= bounds[i]) throw Exception(
		               "Bound " + TextTools::toString(i+1) + " (" + TextTools::toString(bounds[i+1]) +
				") is <= to bound " + TextTools::toString(i)   + " (" + TextTools::toString(bounds[i]) + ").");
    _midPoints[i] = (_bounds[i] + _bounds[i+1]) / 2.;
  }
}

Domain::Domain(Vdouble bounds, Vdouble midPoints) throw (Exception): _bounds(bounds)
{
  if(bounds.size() != midPoints.size() + 1) throw Exception("Domain::Domain(). Number of midpoints must equal number of bounds - 1.");
  
  unsigned int n = _bounds.size();
  for(unsigned int i = 0; i < n-1; i++) {
		if(bounds[i+1] <= bounds[i]) throw Exception(
  	               "Bound " + TextTools::toString(i+1) + " (" + TextTools::toString(bounds[i+1]) +
    	  ") is <= to bound " + TextTools::toString(i)   + " (" + TextTools::toString(bounds[i]) + ").");
  }
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

unsigned int Domain::getIndex(double x) const throw (OutOfRangeException)
{
	if(x < getLowerBound() || x >= getUpperBound())
		throw OutOfRangeException("Domain::getIndex", x, getLowerBound(), getUpperBound());
	for(unsigned int i = 1; i < _bounds.size(); i++)
		if(x < _bounds[i])
			return i - 1;
	// This line can't be reached:
	throw Exception("Unexpected error!");
}

