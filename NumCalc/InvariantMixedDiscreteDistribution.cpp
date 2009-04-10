//
// File: InvariantMixedDiscreteDistribution.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 24 12:02 2007
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

#include "InvariantMixedDiscreteDistribution.h"

//From Utils:
#include <Utils/MapTools.h>

using namespace bpp;

/******************************************************************************/

InvariantMixedDiscreteDistribution::InvariantMixedDiscreteDistribution(
    DiscreteDistribution* dist, double p, double invariant, bool ownDist):
  dist_(dist),
  ownDist_(ownDist),
  invariant_(invariant)
{
  _parameters.addParameters(dist_->getParameters());
  _parameters.addParameter(Parameter("p", p, &Parameter::PROP_CONSTRAINT_IN));
  applyParameters();
}

/******************************************************************************/

void InvariantMixedDiscreteDistribution::fireParameterChanged(const ParameterList & parameters)
{
  dist_->matchParametersValues(parameters);
  applyParameters();
}

/******************************************************************************/


Domain InvariantMixedDiscreteDistribution::getDomain() const
{
  return Domain(bounds_, MapTools::getKeys<double, double, AbstractDiscreteDistribution::Order>(_distribution));
}

/******************************************************************************/

void InvariantMixedDiscreteDistribution::applyParameters()
{
  double p = _parameters.getParameter("p")->getValue();
  _distribution.clear();
  _distribution[invariant_] = p;
  unsigned int distNCat = dist_->getNumberOfCategories();
  bounds_.clear();
  bounds_.push_back(0.);
  vector<double> probs = dist_->getProbabilities();
  vector<double> cats  = dist_->getCategories();
  Domain d = dist_->getDomain();
  for(unsigned int i = 0; i < distNCat; i++)
  {
    _distribution[cats[i]] = (1. - p) * probs[i];
    bounds_.push_back(d.getBound(i));
  }
  bounds_.push_back(d.getBound(distNCat));
}

/******************************************************************************/

