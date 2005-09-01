//
// File RandomTools.cpp
// Author : Julien Dutheil
// Last modification : Friday Septembre 24 2004
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

#include "RandomTools.h"
#include "Uniform01K.h"

// Class destructor
RandomTools::~RandomTools() {}

RandomFactory * RandomTools::DEFAULT_GENERATOR = new Uniform01K(time(NULL));

// Initiate random seed :
//RandomTools::RandInt RandomTools::r = time(NULL) ;

void RandomTools::setSeed(long seed) {
	//r.setSeed(seed);
	DEFAULT_GENERATOR -> setSeed(seed);
}

// Method to get a double random value (between 0 and specified range)
// Note : the number you get is between 0 and entry not including entry !
double RandomTools::giveRandomNumberBetweenZeroAndEntry(double entry, const RandomFactory * generator) {
	//double tm = r.drawFloatNumber();
	double tm = generator -> drawNumber();
	return (tm * entry);
}

// Method to get a boolean random value
bool RandomTools::flipCoin(const RandomFactory * generator) {
	return ((RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator) - 0.5) > 0);
}

// Method to get a integer random value (between 0 and specified range)
// Note : the number you get is between 0 and entry not including entry !
int RandomTools::giveIntRandomNumberBetweenZeroAndEntry(int entry, const RandomFactory * generator) {
	return (int)(giveRandomNumberBetweenZeroAndEntry(entry, generator));
}

double RandomTools::randGaussian(double mean, double variance, const RandomFactory * generator) {
	static int N = 100;
	static double X;
	X=0.0-N/2; /* set mean to 0 */
	for (int ri = 0; ri < N; ri++){
		//    X += 1.0*rand()/RAND_MAX;
		X += giveRandomNumberBetweenZeroAndEntry(1, generator);
	}
	
	/* for uniform randoms in [0,1], mu = 0.5 and var = 1/12 */
	/* adjust X so mu = 0 and var = 1 */
		
	//  X = X * sqrt(12 / N);       /* adjust variance to 1 */
	//  cout <<X * sqrt(variance*12.0/N) + mean<<" ";
	double g = X * sqrt(variance*12.0/N) + mean;
	return (g);
}

double RandomTools::randGamma(double dblAlpha, const RandomFactory * generator) {
	assert(dblAlpha > 0.0);
	if( dblAlpha < 1.0 ) return RandomTools::DblGammaLessThanOne(dblAlpha, generator);
	else if( dblAlpha > 1.0 ) return RandomTools::DblGammaGreaterThanOne(dblAlpha, generator);
	return -log(RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator));
}  

double RandomTools::randGamma(double alpha, double beta, const RandomFactory * generator) {
	double x= RandomTools::randGamma(alpha, generator) / beta;
	return x;
}

double RandomTools::randExponential(double mean, const RandomFactory * generator) {
	return - mean * log(RandomTools::giveRandomNumberBetweenZeroAndEntry(1, generator));
}

//------------------------------------------------------------------------------

	
double RandomTools::DblGammaGreaterThanOne(double dblAlpha, const RandomFactory * generator) {
	// Code adopted from David Heckerman
 	//-----------------------------------------------------------
 	//	DblGammaGreaterThanOne(dblAlpha)
	//
	//	routine to generate a gamma random variable with unit scale and
	//      alpha > 1
 	//	reference: Ripley, Stochastic Simulation, p.90 
 	//	Chang and Feast, Appl.Stat. (28) p.290
 	//-----------------------------------------------------------
    double rgdbl[6];
    
    rgdbl[1] = dblAlpha - 1.0;
    rgdbl[2] = (dblAlpha - (1.0 / (6.0 * dblAlpha))) / rgdbl[1];
    rgdbl[3] = 2.0 / rgdbl[1];
    rgdbl[4] = rgdbl[3] + 2.0;
    rgdbl[5] = 1.0 / sqrt(dblAlpha);
    
    for (;;) {
		double dblRand1;
		double dblRand2;
		do {
	    	dblRand1 = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator);
	    	dblRand2 = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator);
	        if (dblAlpha > 2.5) dblRand1 = dblRand2 + rgdbl[5] * (1.0 - 1.86 * dblRand1);
		} while (!(0.0 < dblRand1 && dblRand1 < 1.0));
	
		double dblTemp = rgdbl[2] * dblRand2 / dblRand1;
	
		if (rgdbl[3] * dblRand1 + dblTemp + 1.0 / dblTemp <= rgdbl[4] ||
	    	rgdbl[3] * log(dblRand1) + dblTemp - log(dblTemp) < 1.0)  {
	    	return dblTemp * rgdbl[1];
		}
    }
    assert(false);
    return 0.0;
}

double RandomTools::DblGammaLessThanOne(double dblAlpha, const RandomFactory * generator) {
	//routine to generate a gamma random variable with 
	//unit scale and alpha < 1
	//reference: Ripley, Stochastic Simulation, p.88 
	double dblTemp;
	const double dblexp = exp(1.0);
	for (;;) {
		double dblRand0 = giveRandomNumberBetweenZeroAndEntry(1.0, generator);
		double dblRand1 = giveRandomNumberBetweenZeroAndEntry(1.0, generator);
		if (dblRand0 <= (dblexp / (dblAlpha + dblexp))){
			dblTemp = pow(((dblAlpha + dblexp) * dblRand0) /
			dblexp, 1.0 / dblAlpha);
			if (dblRand1 <= exp(-1.0 * dblTemp)) return dblTemp;
		} else {
			dblTemp = -1.0 * log((dblAlpha + dblexp) * (1.0 - dblRand0) / (dblAlpha * dblexp)); 
			if (dblRand1 <= pow(dblTemp,dblAlpha - 1.0)) return dblTemp;
		}
	}
    assert(false);
    return 0.0;
}
