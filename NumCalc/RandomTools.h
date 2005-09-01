//
// File RandomTools.h
// Author : Julien Dutheil
//          Sylvain Gaillard
// Last modification : Friday September 24 2004
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

#ifndef _RANDOMTOOLS_H_
#define _RANDOMTOOLS_H_

// From the STL:
#include <cmath>
#include <cassert>
#include <ctime>
#include <vector>
using namespace std;

#include "RandomFactory.h"

class RandomTools
{
	public:
		// Class destructor
		~RandomTools();
	

	public:
		static RandomFactory * DEFAULT_GENERATOR;
		
		/**
		 * @brief Get a double random value (between 0 and specified range).
		 *
		 * Note : the number you get is between 0 and entry not including entry !
		 * @param entry Max number to reach.
		 * @param generator Random number generator to use.
		 */
		static double giveRandomNumberBetweenZeroAndEntry(double entry, const RandomFactory * generator = DEFAULT_GENERATOR);

		/**
		 * @brief Get a boolean random value.
		 *
		 * @param generator Random number generator to use.
		 */
		static bool flipCoin(const RandomFactory * generator = DEFAULT_GENERATOR);

		/**
		 * @brief Get an integer random value (between 0 and specified range).
		 *
		 * Note : the number you get is between 0 and entry not including entry !
		 * @param entry Max number to reach.
		 * @param generator Random number generator to use.
		 */
		static int giveIntRandomNumberBetweenZeroAndEntry(int entry, const RandomFactory * generator = DEFAULT_GENERATOR);

		/**
		 * @brief Set the default generator seed.
		 *
		 * @param seed New seed.
		 */
		static void setSeed(long seed);

		static double randGaussian(double mean, double variance, const RandomFactory * generator = DEFAULT_GENERATOR);
		
		// Routine to generate a gamma random variable with unit scale (beta = 1)
		static double randGamma(double dblAlpha, const RandomFactory * generator = DEFAULT_GENERATOR);

		static double randGamma(double alpha, double beta, const RandomFactory * generator = DEFAULT_GENERATOR);
	
  	static double randExponential(double mean, const RandomFactory * generator = DEFAULT_GENERATOR);

		/**
		 * @brief Get a sample of a vector.
		 *
		 * The sample is a new vector of the specified size.
		 * If the size of the sample is identical to the original vector,
		 * the result is a shuffl of the original vector.
		 *
		 * @param v The vector to sample.
		 * @param size The size of the sample.
		 * @return A vector which is a sample of v.
		 */
		template<class T> static vector<T> getSample(const vector<T> & v, unsigned int size) {
			vector<unsigned int> hat;
			for (unsigned int i = 0 ; i < v.size() ; i++)
				hat.push_back(i);
			vector<T> sample;
			for (unsigned int i = 0 ; i < size ; i++) {
				unsigned int pos = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(hat.size());
				sample.push_back(v[hat[pos]]);
				hat[pos] = hat[hat.size() - 1];
				hat.pop_back();
			}
			return sample;
		}

	private:
		static double DblGammaGreaterThanOne(double dblAlpha, const RandomFactory * generator);
		static double DblGammaLessThanOne(double dblAlpha, const RandomFactory * generator);
};

#endif	//_RANDOMTOOLS_H_

