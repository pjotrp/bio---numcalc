//
// File: ConstantDistribution.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 24 08:48:03 2003
//

#ifndef _CONSTANTDISTRIBUTION_H_
#define _CONSTANTDISTRIBUTION_H_

#include "DiscreteDistribution.h"

class ConstantDistribution : public AbstractDiscreteDistribution
{
    protected:
        double _value;
        
    public:
		ConstantDistribution(double value);
		virtual ~ConstantDistribution() {}
	
	public:
        Domain getDomain() const;
		void fireParameterChanged() {}
	
};


#endif	//_CONSTANTDISTRIBUTION_H_
