/*
 * File RandomTools.h
 * Author : Julien Dutheil <julien.dutheil@ens-lyon.fr>
 * Last modification : Friday September 24 2004
*/

/*
 * This code is adapted from Tal Pupko's SEMPHY library.
*/

// Secured inclusion of header's file
#ifndef _RANDOMTOOLS_H_
#define _RANDOMTOOLS_H_

// From the STL:
#include <cmath>
#include <cassert>
#include <ctime>
#include "RandomFactory.h"

class RandomTools
{
	public:
		// Class destructor
		~RandomTools();
	
// Old Tal's class:
//	public:
//		class RandInt {
//	
//			protected:
//				unsigned long randx;
//
//			public:
//				RandInt(long s = 0) { randx = s; }
//				void setSeed(long s) { randx = s; }
//				int abs(int x) { return x & 0x7fffffff; }
//				static double getMaxNumber() { return 2147483648.0; }
//				int drawNumber() { return randx = randx * 1103515245 + 12345; }
//				double drawFloatNumber() { return abs(drawNumber()) / getMaxNumber(); } //random number between zero and 1
//		};

	public:
		static RandomFactory * DEFAULT_GENERATOR;
		// Method to get a double random value (between 0 and specified range)
		// Note : the number you get is between 0 and entry not including entry !
		static double giveRandomNumberBetweenZeroAndEntry(double entry, const RandomFactory * generator = DEFAULT_GENERATOR);

		// Method to get a boolean random value
		static bool flipCoin(const RandomFactory * generator = DEFAULT_GENERATOR);

		// Method to get a integer random value (between 0 and specified range)
		// Note : the number you get is between 0 and entry not including entry !
		static int giveIntRandomNumberBetweenZeroAndEntry(int entry, const RandomFactory * generator = DEFAULT_GENERATOR);

		static void setSeed(long seed);

		static double randGaussian(double mean, double variance, const RandomFactory * generator = DEFAULT_GENERATOR);
		
		// Routine to generate a gamma random variable with unit scale (beta = 1)
		static double randGamma(double dblAlpha, const RandomFactory * generator = DEFAULT_GENERATOR);

		static double randGamma(double alpha, double beta, const RandomFactory * generator = DEFAULT_GENERATOR);
	
		// Added by jdutheil on 08 11 2002.
  		static double randExponential(double mean, const RandomFactory * generator = DEFAULT_GENERATOR);

	private:
//		static RandInt r;
		static double DblGammaGreaterThanOne(double dblAlpha, const RandomFactory * generator);
		static double DblGammaLessThanOne(double dblAlpha, const RandomFactory * generator);
};

#endif	//_RANDOMTOOLS_H_

