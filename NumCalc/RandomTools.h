/*
 * File RandomTools.h
 * Author : Julien Dutheil <julien.dutheil@ens-lyon.fr>
 * Last modification : Tuesday August 21 2003
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

/* Static class RandomTools
 * Implements methods to get random numbers from different algorithms
 * 
 * Initial value for randomize (random seed) is managed in methods
 *
 * This class can't be instanciated
*/

class RandomTools
{
	public:
		// Class destructor
		~RandomTools();
	
	public:
		class RandInt {
	
			protected:
				unsigned long randx;

			public:
				RandInt(long s = 0) { randx = s; }
				void setSeed(long s) { randx = s; }
				int abs(int x) { return x & 0x7fffffff; }
				static double getMaxNumber() { return 2147483648.0; }
				int drawNumber() { return randx = randx * 1103515245 + 12345; }
				double drawFloatNumber() { return abs(drawNumber()) / getMaxNumber(); } //random number between zero and 1
		};

	public:
		// Method to get a double random value (between 0 and specified range)
		// Note : the number you get is between 0 and entry not including entry !
		static double giveRandomNumberBetweenZeroAndEntry(double entry);

		// Method to get a boolean random value
		static bool flipCoin();

		// Method to get a integer random value (between 0 and specified range)
		// Note : the number you get is between 0 and entry not including entry !
		static int giveIntRandomNumberBetweenZeroAndEntry(int entry);

		static void setSeed(unsigned long seed);

		static double randGaussian(double mean, double variance);
		
		// Routine to generate a gamma random variable with unit scale (beta = 1)
		static double randGamma(double dblAlpha);

		static double randGamma(double alpha, double beta);
	
		// Added by jdutheil on 08 11 2002.
  		static double randExponential(double mean);

	private:
		static RandInt r;
		static double DblGammaGreaterThanOne(double dblAlpha);
		static double DblGammaLessThanOne(double dblAlpha);
};

#endif	//_RANDOMTOOLS_H_
