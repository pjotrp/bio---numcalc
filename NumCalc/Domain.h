//
// File: Domain.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Feb  4 19:12:44 2004
//

#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "VectorTools.h"

// From Utils:
#include <Utils/Exceptions.h>
#include <Utils/TextTools.h>

class OutOfRangeException: public Exception
{
	protected:
		double _badValue, _lowerBound, _upperBound;
	
	public:
		OutOfRangeException(const string & text, double badValue, double lowerBound, double upperBound):
			Exception(text +
				"Value " +
				TextTools::toString(badValue) + 
				" is out of range [" +
				TextTools::toString(lowerBound) + 
				", " + 
				TextTools::toString(upperBound) +
				"[."),
			_badValue(badValue), _lowerBound(lowerBound), _upperBound(upperBound) {}
};


/**
 * @brief This class is quite a C++ implementation of the Java PSOL library.
 */ 
class Domain
{
	protected:
		Vdouble _bounds;
		Vdouble _midPoints;
	
	public:
		Domain(double a, double b, unsigned int n);
        Domain(Vdouble bounds);
        Domain(Vdouble bounds, Vdouble midPoints) throw (Exception);
	
		virtual ~Domain() {}
	
	public:
		virtual double getLowerBound() const;
		virtual double getUpperBound() const;
		virtual double getLowerValue() const;
		virtual double getUpperValue() const;
		virtual unsigned int getSize() const;
		virtual double getBound(unsigned int i) const;
		virtual double getValue(unsigned int i) const;
		virtual double getNearestValue(double x) const throw (OutOfRangeException);
		virtual int getIndex(double x) const throw (OutOfRangeException);
		
};


#endif	//_DOMAIN_H_
