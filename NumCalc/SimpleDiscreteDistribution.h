#ifndef _SIMPLEDISCRETEDISTRIBUTION_H_
#define _SIMPLEDISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"

// From the STL:
#include <map>
using namespace std;

class SimpleDiscreteDistribution: public AbstractDiscreteDistribution {
	public:
		SimpleDiscreteDistribution(
				const map<double, double> & distribution);
		virtual ~SimpleDiscreteDistribution();

	public:
		void fireParameterChanged() {}
		Domain getDomain() const;
};

#endif	//_SIMPLEDISCRETEDISTRIBUTION_H_
