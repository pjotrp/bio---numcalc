//
// File: DiscreteDistribution.h
// Created by:  <>
// Created on: Sat Apr 19 14:37:30 2003
//

#ifndef _DISCRETEDISTRIBUTION_H_
#define _DISCRETEDISTRIBUTION_H_

#include <map>

using namespace std;

#include "VectorTools.h"
#include "Parametrizable.h"
#include "Domain.h"

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

};

class AbstractDiscreteDistribution: public DiscreteDistribution {

    public:
    
        class Order {
            public :
                bool operator() (double l1, double l2) const {
                    return (l1 < l2 - 0.000000001); //precision of E9
                }
        };

	protected:
		//these fields must be initialized in the constructor of the derived classes:
		map<double, double, Order> _distribution;
		mutable ParameterList _parameters;
	
	public:
		AbstractDiscreteDistribution() {}
		virtual ~AbstractDiscreteDistribution() {}
	
	public:
        unsigned int getNumberOfCategories() const;
		double getCategory(unsigned int categoryIndex) const;
		double getProbability(unsigned int categoryIndex) const;
		double getProbability(double category) const;
		Vdouble getCategories() const;
		Vdouble getProbabilities() const;
		void set(double category, double probability);
		void add(double category, double probability);
		//Pr(x < category):
		double  getInfCumulativeProbability(double category) const;
		//Pr(x <= category):
		double getIInfCumulativeProbability(double category) const;
		//Pr(x > category):
		double  getSupCumulativeProbability(double category) const;
		//Pr(x >= category):
		double getSSupCumulativeProbability(double category) const;

		double rand() const;

		virtual void fireParameterChanged() = 0;
		
		//Parametrizable interface implemented here:
		ParameterList getParameters()                            const ;
		double        getParameter (const string & name)         const throw (ParameterNotFoundException);
		void       setAllParametersValues(const ParameterList & params)      throw (ParameterNotFoundException, ConstraintException);
		void          setParameterValue  (const string & name, double value) throw (ParameterNotFoundException, ConstraintException);
		void          setParametersValues(const ParameterList & params)      throw (ParameterNotFoundException, ConstraintException);
		void        matchParametersValues(const ParameterList & params)      throw (ConstraintException);
};

#endif	//_DISCRETEDISTRIBUTION_H_
