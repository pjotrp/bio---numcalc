#ifndef _DISCRETEDISTRIBUTION_H_
#define _DISCRETEDISTRIBUTION_H_

#include "VectorTools.h"
#include "Parametrizable.h"
#include "Domain.h"

// From the STL:
#include <iostream>
using namespace std;

class DiscreteDistribution: public Parametrizable {
	public:
		DiscreteDistribution() {}
		virtual ~DiscreteDistribution() {}
	
	public:
		virtual unsigned int getNumberOfCategories() const = 0;
		virtual double getCategory(unsigned int categoryIndex) const = 0;
		virtual double getProbability(unsigned int categoryIndex) const = 0;
		virtual double getProbability(double category) const = 0;
		virtual Vdouble getCategories() const = 0;
		virtual Vdouble getProbabilities() const = 0;
		virtual void set(double category, double probability) = 0;
		virtual void add(double category, double probability) = 0;
		//Pr(x < category):
		virtual double  getInfCumulativeProbability(double category) const = 0;
		//Pr(x <= category):
		virtual double getIInfCumulativeProbability(double category) const = 0;
		//Pr(x > category):
		virtual double  getSupCumulativeProbability(double category) const = 0;
		//Pr(x >= category):
		virtual double getSSupCumulativeProbability(double category) const = 0;
	
		virtual double rand() const = 0;

    virtual Domain getDomain() const = 0;

		/**
		 * @brief Print the distribution (categories and corresponding probabilities).
		 *
		 * @param out The outstream where to print the distribution.
		 */
		virtual void print(ostream & out) const = 0;

};

#endif	//_DISCRETEDISTRIBUTION_H_
