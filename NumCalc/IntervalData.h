//
// File: IntervalData.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Feb  4 19:59:29 2004
//

#ifndef _INTERVALDATA_H_
#define _INTERVALDATA_H_

//#include "VectorTools.h"
#include "Domain.h"

// from the STL:
#include <iostream>
using namespace std;

/**
 * @brief This class is quite a C++ implementation of the Java PSOL library.
 */ 
class IntervalData
{
	protected:
		//Vdouble _data;
		const Domain _domain;
		vector<unsigned int> _freqs;
		string _name;
		double _sum, _sumsquare, _min, _max;
		unsigned int _n;
	
	public:
		IntervalData(const Domain & domain, const string & name = "");
	
		virtual ~IntervalData() {}
			
	public:
		virtual Domain getDomain() const;
		virtual double getDomainValue(double x) const;
		virtual void setName(const string & name);
		virtual string getName() const;
		virtual unsigned int getFreq(unsigned int i) const;
		virtual double getDensity(unsigned int i) const;
		virtual vector<unsigned int> getFrequencies() const;
		virtual Vdouble getDensities() const;
		virtual void addValue(double value) throw (OutOfRangeException);
		virtual unsigned int getSize() const;
		virtual double getMinValue() const;
		virtual double getMaxValue() const;
	
		/**
		 * @brief Return the true mean of the dataset.
		 */
		virtual double getMean() const;
		
		/**
		 * @brief Return the standard deviation of the data set treated as a
		 * sample (that is, the sum is divided by @f$ n - 1 @f$ where @f$ n @f$
		 * is the number of points.
		 */
		virtual double getSD() const;
	
		/**
		 * @brief Return the standard deviation of the data set treated as a
		 * population (that is, the sum is divided by the number of points @f$n @f$
		 * rather than @f$ n-1 @f$ where @f$ n @f$  is the number of points.
		 */
		virtual double getSDP() const;
		
		virtual void reset();
		
		virtual void print(ostream & out) const;

};


#endif	//_INTERVALDATA_H_
