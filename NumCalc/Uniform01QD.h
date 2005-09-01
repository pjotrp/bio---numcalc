//
// File Uniform01QD.h
// Author : Sylvain Gaillard
// Last modification : Friday September 24 2004
//
//
/*
Copyright or � or Copr. CNRS, (November 17, 2004)

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


#ifndef _UNIFORM01QD_H_
#define _UNIFORM01QD_H_

#include "RandomFactory.h"

/**
 * @brief A quick and dirty uniform random number generator.
 *
 * This is a congruential uniform generator which draw double between 0 and 1 excluding
 * the end points.
 *
 * WARNING!!! Only works on 32bits architectures!
 */
class Uniform01QD : public RandomFactory {
	public: // Constructors and destructor
		/**
		 * @brief Create a Random Number Generator.
		 *
		 * @param seed The seed for the random numbers.
		 */
		Uniform01QD(long seed);

		/**
		 * @brief Destroy the generator.
		 */
		~Uniform01QD();

	public:
		/**
		 * @brief Set the seed for a new set of random numbers.
		 */
		void setSeed(long seed);

		/**
		 * @brief Get a random number between 0.0 and 1.0 (exclusive of the end point values).
		 */
		double drawNumber() const;

	protected:
		mutable long _seed;
};
#endif // _UNIFORM01QD_H_
