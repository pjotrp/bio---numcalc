//
// File: Constraints.h
// Created by: Julien Dutheil
// Created on: Thu Dec 25 19:35:17 2003
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

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

#ifndef _CONSTRAINTS_H_
#define _CONSTRAINTS_H_

// From the STL:
#include <string>
using namespace std;

//From Utils:
#include <Utils/TextTools.h>

/**
 * @brief The constraint interface.
 *
 * It provides a method that tells if a given value is correct.
 */
class Constraint {
	
	public:
		Constraint() {}
		virtual ~Constraint() {}
			
	public:
		/**
		 * @brief Tell if a given value is correct.
		 *
		 * @param value The value to test.
		 * @return True is the value is correct.
		 */
		virtual bool isCorrect(double value) const = 0;
			
		/**
		 * @brief Give the nearest limit for a bad value.
		 *
		 * @param value The bad value.
		 * @return The nearer limit.
		 */
		virtual double getLimit(double value) const = 0;

		/**
		 * @brief Give a short description on the type of constraint.
		 *
		 * @return A string which describes the constraint.
		 */
		virtual string getDescription() const = 0;		
};
	
/**
 * @brief Including positive real constraint.
 */
class IncludingPositiveReal: public Constraint
{
	protected:
		double _lower;
	
	public:
		IncludingPositiveReal(double lowerBound): _lower(lowerBound) {}
		virtual ~IncludingPositiveReal() {}
	
	public:
		bool isCorrect(double value) const { return value >= _lower; }
		double getLimit(double value) const { return isCorrect(value) ? value : _lower; }
		string getDescription() const
		{
			return "[ " + TextTools::toString(_lower) + ", +inf [";
		}
		
};
		
/**
 * @brief Excluding positive real constraint.
 */
class ExcludingPositiveReal: public Constraint
{
	protected:
		double _lower;
	
	public:
		ExcludingPositiveReal(double lowerBound): _lower(lowerBound) {}
		virtual ~ExcludingPositiveReal() {}
	
	public:
		bool isCorrect(double value) const { return value > _lower; }
		double getLimit(double value) const { return isCorrect(value) ? value : _lower; }
		string getDescription() const
		{
			return "] " + TextTools::toString(_lower) + ", +inf [";
		}
};

/**
 * @brief Interval constraint base class.
 */
class Interval: public Constraint
{
	protected:
		double _lower, _upper;
		
	public:
		/**
		 * @brief Build a new interval constraint.
		 *
		 * @param lowerBound The lower bound of the interval.
		 * @param upperBound The upper bound of the interval.
		 */
		Interval(double lowerBound, double upperBound): _lower(lowerBound), _upper(upperBound) {}
		
		virtual ~Interval() {}
			
	public:
		bool isCorrect(double value) const = 0;
		double getLimit(double value) const
		{
			if(isCorrect(value)) return value;
			else if (value <= _lower) return _lower;
			else return _upper;
		}	
		string getDescription() const = 0;
};

/**
 * @brief Including interval.
 */
class IncludingInterval: public Interval
{
	public:
		/**
		 * @brief Build a new including interval constraint.
		 *
		 * @param lowerBound The lower bound of the interval.
		 * @param upperBound The upper bound of the interval.
		 */
		IncludingInterval(double lowerBound, double upperBound) :
			Interval(lowerBound, upperBound) {}

		virtual ~IncludingInterval() {}
	
	public:
		bool isCorrect(double value) const
		{
			return value >= _lower && value <= _upper;
		}
		string getDescription() const
		{
			return "[ " + TextTools::toString(_lower) + ", "
	  		          + TextTools::toString(_upper) + " ]";
		}
};
		
/**
 * @brief Excluding interval.
 */
class ExcludingInterval: public Interval
{
	public:
		/**
		 * @brief Build a new excluding interval constraint.
		 *
		 * @param lowerBound The lower bound of the interval.
		 * @param upperBound The upper bound of the interval.
		 */
		ExcludingInterval(double lowerBound, double upperBound):
			Interval(lowerBound, upperBound) {}
    virtual ~ExcludingInterval() {}
	
	public:
		bool isCorrect(double value) const
		{
			return value > _lower && value < _upper;
		}
		string getDescription() const
		{
			return "] " + TextTools::toString(_lower) + ", "
    		         	+ TextTools::toString(_upper) + " [";
		}
		
};
		
/**
 * @brief Left-including, right-excluding interval.
 */
class IncludingExcludingInterval: public Interval
{
	public:
		/**
		 * @brief Build a new left-including, right-excluding interval constraint.
		 *
		 * @param lowerBound The lower bound of the interval.
		 * @param upperBound The upper bound of the interval.
		 */
		IncludingExcludingInterval(double lowerBound, double upperBound):
			Interval(lowerBound, upperBound) {}
		
    virtual ~IncludingExcludingInterval() {}
	
	public:
		bool isCorrect(double value) const
		{
			return value >= _lower && value < _upper;
		}

		string getDescription() const
		{
			return "[ " + TextTools::toString(_lower) + ", " 
	  		          + TextTools::toString(_upper) + " [";
		}
		
};
		
/**
 * @brief Left-excluding, right-including interval.
 */
class ExcludingIncludingInterval: public Interval
{
	public:
		/**
		 * @brief Build a new left-excluding, right-including interval constraint.
		 *
		 * @param lowerBound The lower bound of the interval.
		 * @param upperBound The upper bound of the interval.
		 */
		ExcludingIncludingInterval(double lowerBound, double upperBound) :
			Interval(lowerBound, upperBound) {}

		virtual ~ExcludingIncludingInterval() {}

	public:
		bool isCorrect(double value) const
		{
			return value > _lower && value <= _upper;
		}
		string getDescription() const
		{
			return "] " + TextTools::toString(_lower) + ", "
      		      	+ TextTools::toString(_upper) + " ]";
		}
};

#endif	//_CONSTRAINTS_H_

