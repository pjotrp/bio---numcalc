/*
 * File Uniform01QD.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
 */

// Secured inclusion of header's file
#ifndef _UNIFORM01QD_H_
#define _UNIFORM01QD_H_

#include "RandomFactory.h"

/**
 * @brief A quick and dirty uniform random number generator.
 *
 * This is a congruential uniform generator which draw double between 0 and 1 excluding
 * the end points.
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
