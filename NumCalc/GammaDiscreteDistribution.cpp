//
// File: GammaDiscreteDistribution.cpp
// Created by: Julien Dutheil
// Created on: Sun Oct 26 20:36:12 2003
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

#include "GammaDiscreteDistribution.h"
#include "RandomTools.h"

// From Utils:
#include <Utils/MapTools.h>

// From the STL:
#include <cmath>

const double GammaDiscreteDistribution::VERYBIG = static_cast<double>(1.7E+23);

/** Constructor: **************************************************************/

GammaDiscreteDistribution::GammaDiscreteDistribution(unsigned int n, double alpha) : AbstractDiscreteDistribution()
{
	// We use a lower bound of 0.05 for alpha to prohibe errors due to computer
	// floating precision: if alpha is quite low (gamma -> constant), some classes
	// may have the same category value, leading to a classe number lower than expected.
	// NB: if this is the case, then a warning is shown. This may happen in optimization
	// algorithms.
	_alphaConstraint = new IncludingPositiveReal(0.05);
	_parameters.addParameter(Parameter("alpha", alpha, _alphaConstraint));
	applyParameters(n);
}


GammaDiscreteDistribution::~GammaDiscreteDistribution() {
	delete _alphaConstraint;
}

/******************************************************************************/

void GammaDiscreteDistribution::fireParameterChanged(const ParameterList & parameters)
{
	//cout << "Parameter changed, alpha = " << getParameter("alpha") << endl;
	applyParameters(getNumberOfCategories());	
}

/******************************************************************************/

Domain GammaDiscreteDistribution::getDomain() const
{
    return Domain(_bounds, MapTools::getKeys<double, double, AbstractDiscreteDistribution::Order>(_distribution));
}

/******************************************************************************/

void GammaDiscreteDistribution::applyParameters(unsigned int numberOfCategories) {
	if(numberOfCategories <= 0)
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in GammaDiscreteDistribution::applyParameters()." << endl;
	_distribution.clear();
	_bounds.clear();
	_bounds.resize(numberOfCategories + 1);
	//cout<<"number of categories is: "<<categories()<<endl;
	//cout<<"alpha is: "<<alpha<<endl;
	if (numberOfCategories == 1) {
		_distribution[0] = 1.0;
		_bounds[0] = 0; _bounds[1] = VERYBIG;
		return;
	} else if (numberOfCategories > 1) {
		double alpha = getParameterValue("alpha");
		vector<double> means = computeValues(numberOfCategories, alpha, alpha, false);
		double p = 1. / (double)numberOfCategories;
		for(unsigned int i = 0; i < numberOfCategories; i++) {
			_distribution[means[i]] = p;
		}
		_bounds = computeBounds(numberOfCategories, alpha, alpha);
		if(getNumberOfCategories() != numberOfCategories) {
			cout << "WARNING!!! Couldn't create " << numberOfCategories << " distinct categories." << endl;
			cout << "WARNING!!! This may occure if you specified a too low alpha parameter." << endl;
		}
		return ;
	} else {
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in GammaDiscreteDistribution::applyParameters()." << endl;
	}
}

/******************************************************************************/

// Adapted from function DiscreteGamma of Yang

vector<double> GammaDiscreteDistribution::computeValues(unsigned int nbClasses, double alfa, double beta, bool median)
{
/* discretization of gamma distribution with equal proportions in each 
   category
*/
   unsigned int K = nbClasses;
   vector<double> rK(K), freqK(K);
	 int i;
	 double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

   if (median) {
      for (i=0; i<K; i++) rK[i]=RandomTools::qGamma((i*2.0+1)*gap05, alfa, beta);
      for (i=0,t=0; i<K; i++) t+=rK[i];
      for (i=0; i<K; i++)     rK[i]*=factor/t;
   }
   else {
      lnga1=RandomTools::lnGamma(alfa+1);
      for (i=0; i<K-1; i++)
	 freqK[i]=RandomTools::qGamma((i+1.0)/K, alfa, beta);
      for (i=0; i<K-1; i++)
	 freqK[i]=RandomTools::incompleteGamma(freqK[i]*beta, alfa+1, lnga1);
      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return rK;
}

vector<double> GammaDiscreteDistribution::computeBounds(unsigned int nbClasses, double alfa, double beta)
{
	vector<double> bounds(nbClasses + 1);
	bounds[0] = 0;
	for(unsigned int i = 1; i < nbClasses; i++) bounds[i] = RandomTools::qGamma((double)i/(double)nbClasses, alfa, beta);
	bounds[nbClasses] = VERYBIG;
	return bounds;
}

