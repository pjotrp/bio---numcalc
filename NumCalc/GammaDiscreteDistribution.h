//
// File: GammaDiscreteDistribution.h
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

#ifndef _GAMMADISCRETEDISTRIBUTION_H_
#define _GAMMADISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"
#include "Constraints.h"
#include "RandomTools.h"

/**
 * @brief Discretized Gamma distribution.
 */
class GammaDiscreteDistribution:
  public AbstractDiscreteDistribution
{
	protected:
		vector<double> _bounds;
		IncludingPositiveReal * _alphaConstraint;
	
		static const double VERYBIG;
		
	public:
		GammaDiscreteDistribution(unsigned int n, double alpha = 1.);

    GammaDiscreteDistribution(const GammaDiscreteDistribution & dist):
      AbstractParametrizable(dist),
      AbstractDiscreteDistribution(dist),
      _bounds(dist._bounds),
      _alphaConstraint(dist._alphaConstraint->clone()) {}
    
    GammaDiscreteDistribution & operator=(const GammaDiscreteDistribution & dist)
    {
      AbstractDiscreteDistribution::operator=(dist);
      _bounds = dist._bounds;
      _alphaConstraint = dist._alphaConstraint->clone();
      return *this;
    }

		virtual ~GammaDiscreteDistribution();

#if defined(VIRTUAL_COV)
    GammaDiscreteDistribution * clone() const
    {
      return new GammaDiscreteDistribution(*this);
    }
#else
    Clonable * clone() const
    {
      return new GammaDiscreteDistribution(*this);
    }
#endif
	
	public:
    Domain getDomain() const;
		void fireParameterChanged(const ParameterList & parameters);
	
  public:
    double randC() const throw (Exception)
    {
      return RandomTools::randGamma(_parameters.getParameter("alpha")->getValue());
    }
    
	protected:
		void applyParameters(unsigned int numberOfCategories);
		static vector<double> computeBounds(unsigned int nbClasses, double alfa, double beta);
		static vector<double> computeValues(unsigned int nbClasses, double alfa, double beta, bool median);
		
};

#endif	//_GAMMADISCRETEDISTRIBUTION_H_

