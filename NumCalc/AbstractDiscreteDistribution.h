#ifndef _ABSTRACTDISCRETEDISTRIBUTION_H_
#define _ABSTRACTDISCRETEDISTRIBUTION_H_

#include "DiscreteDistribution.h"

#include <map>

using namespace std;

class AbstractDiscreteDistribution: public DiscreteDistribution {

  public:
    
    class Order {
      public:
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

		/**
		 * @name From the DiscreteDistribution inteface.
		 *
		 * @{
		 */
    unsigned int getNumberOfCategories() const;
		double getCategory(unsigned int categoryIndex) const;
		double getProbability(unsigned int categoryIndex) const;
		double getProbability(double category) const;
		Vdouble getCategories() const;
		Vdouble getProbabilities() const;
		void set(double category, double probability);
		void add(double category, double probability);
		double  getInfCumulativeProbability(double category) const;
		double getIInfCumulativeProbability(double category) const;
		double  getSupCumulativeProbability(double category) const;
		double getSSupCumulativeProbability(double category) const;
		double rand() const;
		void print(ostream & out) const;
		/** @} */
		
		virtual void fireParameterChanged() = 0;
		
		/**
		 * @name Parametrizable interface implementation.
		 *
		 * @{
		 */
		ParameterList getParameters()                            const ;
		double        getParameter (const string & name)         const throw (ParameterNotFoundException);
		void       setAllParametersValues(const ParameterList & params)      throw (ParameterNotFoundException, ConstraintException);
		void          setParameterValue  (const string & name, double value) throw (ParameterNotFoundException, ConstraintException);
		void          setParametersValues(const ParameterList & params)      throw (ParameterNotFoundException, ConstraintException);
		void        matchParametersValues(const ParameterList & params)      throw (ConstraintException);
		/** @} */
};

#endif	//_ABSTRACTDISCRETEDISTRIBUTION_H_
