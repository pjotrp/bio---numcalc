//
// File: GammaDiscreteDistribution.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Oct 26 20:36:12 2003
//

#ifndef _GAMMADISCRETEDISTRIBUTION_H_
#define _GAMMADISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"
#include "Constraints.h"

#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

class GammaDiscreteDistribution : public AbstractDiscreteDistribution
{
	protected:
		//double _alpha;
		vector<double> _bounds;
		IncludingPositiveReal * _alphaConstraint;
	
		static const double VERYBIG;
		
	public:
		GammaDiscreteDistribution(unsigned int n, double alpha = 1.);
		virtual ~GammaDiscreteDistribution();
	
	public:
    Domain getDomain() const;
		void fireParameterChanged();
	
	protected:
	
		void applyParameters(unsigned int numberOfCategories);
		static vector<double> computeBounds(unsigned int nbClasses, double alfa, double beta);
		static vector<double> computeValues(unsigned int nbClasses, double alfa, double beta, bool median);
		
	//From Yang's PAML package:
		static double PointNormal (double prob);
		static double LnGamma (double alpha);
		static double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
		static double PointChi2(double prob, double v);
};


#endif	//_GAMMADISCRETEDISTRIBUTION_H_
