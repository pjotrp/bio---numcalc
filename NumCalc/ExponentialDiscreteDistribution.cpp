//
// File: ExponentialDiscreteDistribution.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 13 12:37 2007
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "ExponentialDiscreteDistribution.h"
#include "RandomTools.h"

// From Utils:
#include <Utils/MapTools.h>
using namespace bpp;

// From the STL:
#include <cmath>
using namespace std;
  
const double ExponentialDiscreteDistribution::VERYBIG = static_cast<double>(1.7E+23);

/** Constructor: **************************************************************/

ExponentialDiscreteDistribution::ExponentialDiscreteDistribution(unsigned int n, double lambda) : AbstractDiscreteDistribution()
{
	_lambdaConstraint = new IncludingPositiveReal(0.0);
	_parameters.addParameter(Parameter("lambda", lambda, _lambdaConstraint, true));
	applyParameters(n);
}

ExponentialDiscreteDistribution::~ExponentialDiscreteDistribution() {}

/******************************************************************************/

void ExponentialDiscreteDistribution::fireParameterChanged(const ParameterList & parameters)
{
	//cout << "Parameter changed, alpha = " << getParameter("alpha") << endl;
	applyParameters(getNumberOfCategories());	
}

/******************************************************************************/

Domain ExponentialDiscreteDistribution::getDomain() const
{
    return Domain(_bounds, MapTools::getKeys<double, double, AbstractDiscreteDistribution::Order>(_distribution));
}

/******************************************************************************/

void ExponentialDiscreteDistribution::applyParameters(unsigned int numberOfCategories)
{
  discretize(numberOfCategories, _parameters[0]->getValue(), true);
}

/******************************************************************************/

void ExponentialDiscreteDistribution::discretize(unsigned int numberOfCategories, double lambda, bool median)
{
	if(numberOfCategories == 0)
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in ExponentialDiscreteDistribution::applyParameters()." << endl;
	_distribution.clear();
	_bounds.resize(numberOfCategories + 1);
	if(numberOfCategories == 1)
  {
		_distribution[0] = 1.0;
		_bounds[0] = 0; _bounds[1] = VERYBIG;
		return;
	}
  else if(numberOfCategories > 1)
  {
    _bounds.resize(numberOfCategories + 1);
    _distribution.clear();

	  _bounds[0] = 0;
    vector<double> values(numberOfCategories);

	  for(unsigned int i = 1; i <= numberOfCategories; i++)
    {
      double a = _bounds[i-1];
      double b = (i == numberOfCategories)
        ? VERYBIG
        : (1. / lambda) * log((double)numberOfCategories / ((double)(numberOfCategories - i)));
      _bounds[i] = b;
      if(median)
        values[i-1] = (1. / lambda) * log((double)(2*numberOfCategories) / (double)(2*(numberOfCategories - i) + 1)); 
      else
        values[i-1] = (a + 1./lambda) * exp(-a * lambda) + (b + 1. / lambda) * exp(-b * lambda); 
    }

		double p = 1. / (double)numberOfCategories;
		for(unsigned int i = 0; i < numberOfCategories; i++)
    {
			_distribution[values[i]] += p;
		}
		if(getNumberOfCategories() != numberOfCategories)
    {
			cout << "WARNING!!! Couldn't create " << numberOfCategories << " distinct categories." << endl;
			cout << "WARNING!!! This may occure if you specified a too low lambda parameter." << endl;
		}
		return ;
	}
  else
  {
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in ExponentialDiscreteDistribution::applyParameters()." << endl;
	}
}

