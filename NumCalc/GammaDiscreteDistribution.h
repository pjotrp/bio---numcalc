//
// File: GammaDiscreteDistribution.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Oct 26 20:36:12 2003
//

#ifndef _GAMMADISCRETEDISTRIBUTION_H_
#define _GAMMADISCRETEDISTRIBUTION_H_

#include "DiscreteDistribution.h"
#include "Constraints.h"

class GammaDiscreteDistribution : public AbstractDiscreteDistribution
{
	protected:
		//double _alpha;
		vector<double> _bonderi;
		IncludingPositiveReal * _alphaConstraint;
	
		static const double VERYBIG;
		static const int    ITMAX;
		static const double EPS;
		static const double FPMIN;
		static const double ERR_FOR_GAMMA_CALC;
		
	public:
		GammaDiscreteDistribution(unsigned int n, double alpha = 1.);
		virtual ~GammaDiscreteDistribution();
	
	public:
        Domain getDomain() const;
		void fireParameterChanged();
	
	protected:
		//From Tal's Semphy library:
		void   applyParameters(unsigned int numberOfCategories);
		void       fillBonderi(unsigned int numberOfCategories);
		double gammp(double a, double x) throw (BadNumberException);
		double gammln(double xx);

		double search_for_z_in_dis_with_beta_1(double alpha, double ahoson) throw (BadNumberException);
		double search_for_z_in_dis_with_any_beta(double alpha, double beta, double ahoson);
		void   gser(double *gamser, double a, double x, double *gln) throw (BadNumberException);
		void   gcf (double *gammcf, double a, double x, double *gln) throw (BadNumberException);
		double the_average_r_in_category_between_a_and_b(double a, double b, double alpha, double beta, int k);

		//This method was the one used in the Semphy library:
		void fillMean(unsigned int numberOfCategories);
		//But I replaced it by this one, using Yang's code, in order to follow
		//Yang and use the median of each class and not the mean:
		void fillMedian(unsigned int numberOfCategories);

		//From Yang's PAML package:
		double PointNormal (double prob);
		double CDFNormal (double x);
		double LnGamma (double alpha) throw (BadNumberException);
		double DFGamma(double x, double alpha, double beta) throw (BadNumberException);
		double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
		double PointChi2 (double prob, double v) throw (Exception);
		double PointGamma(double prob, double alpha, double beta) throw (Exception);		
};


#endif	//_GAMMADISCRETEDISTRIBUTION_H_
