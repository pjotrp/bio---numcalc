/*
 * File Uniform01WH.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
 */

// Secured inclusion of header's file
#ifndef _UNIFORM01WH_H_
#define _UNIFORM01WH_H_

#include "RandomFactory.h"

/**
 * @brief A uniform random number generator.
 *
 * This is a congruential uniform generator which draw double between 0 and 1 excluding
 * the end points.
 * This generator is based on a Fortan routine from Wichmann, B. A. and Hill, I. D. (1982).
 * "An efficient and portable pseudorandom number generator," Applied Statistics, 31, 188-190
 */
class Uniform01WH : public RandomFactory {
	public: // Constructors and destructor
		/**
		 * @brief Create a Random Number Generator.
		 *
		 * @param seed The seed for the random numbers.
		 */
		Uniform01WH(long seed);

		/**
		 * @brief Destroy the generator.
		 */
		~Uniform01WH();

	public:
		/**
		 * @brief Set the seed for a new set of random numbers.
		 */
		void setSeed(long seed);

		/**
		 * @brief Set the three seeds.
		 */
		void setSeeds(long seed1, long seed2 = 20356, long seed3 = 35412);

		/**
		 * @brief Get a random number between 0.0 and 1.0 (exclusive of the end point values).
		 */
		double drawNumber() const;

	protected:
		mutable long ix, iy, iz;
};
#endif // _UNIFORM01WH_H_
