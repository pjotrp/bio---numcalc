//
// File: ConstantDistribution.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 24 08:48:03 2003
//

#include "ConstantDistribution.h"

/******************************************************************************/

ConstantDistribution::ConstantDistribution(double value) : AbstractDiscreteDistribution(), _value(value)
{
    _distribution[_value] = 1; //One single class of rate 1 with probability 1.
}

/******************************************************************************/

Domain ConstantDistribution::getDomain() const
{
     return Domain(_value, _value, 1);
}

/******************************************************************************/
