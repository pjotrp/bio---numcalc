//
// File: Constraints.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Dec 25 19:35:17 2003
//

#ifndef _CONSTRAINTS_H_
#define _CONSTRAINTS_H_

// From the STL:
#include <string>
using namespace std;

/**
 * @brief The constraint interface.
 *
 * It provides a single method that tells if a certain value is correct.
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
		 * @brief Give the nearer limit for a bad value.
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
		IncludingPositiveReal(double lowerBound);
		virtual ~IncludingPositiveReal() {}
	
	public:
		bool isCorrect(double value) const;
		double getLimit(double value) const;
		string getDescription() const;
};
		
/**
 * @brief Excluding positive real constraint.
 */
class ExcludingPositiveReal: public Constraint
{
	protected:
		double _lower;
	
	public:
		ExcludingPositiveReal(double lowerBound);
		virtual ~ExcludingPositiveReal() {}
	
	public:
		bool isCorrect(double value) const;
		double getLimit(double value) const;
		string getDescription() const;
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
		Interval(double lowerBound, double upperBound);
		
		virtual ~Interval() {}
			
	public:
		bool isCorrect(double value) const = 0;
		double getLimit(double value) const;
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
		IncludingInterval(double lowerBound, double upperBound);

		virtual ~IncludingInterval() {}
	
	public:
		bool isCorrect(double value) const;
		string getDescription() const;
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
		ExcludingInterval(double lowerBound, double upperBound);
    virtual ~ExcludingInterval() {}
	
	public:
		bool isCorrect(double value) const;
		string getDescription() const;
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
		IncludingExcludingInterval(double lowerBound, double upperBound);
    virtual ~IncludingExcludingInterval() {}
	
	public:
		bool isCorrect(double value) const;
		string getDescription() const;
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
		ExcludingIncludingInterval(double lowerBound, double upperBound);
		virtual ~ExcludingIncludingInterval() {}

	public:
		bool isCorrect(double value) const;
		string getDescription() const;
};


#endif	//_CONSTRAINTS_H_
