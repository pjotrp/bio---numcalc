/*
 * File Uniform01K.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
 */

// Secured inclusion of header's file
#ifndef _UNIFORM01K_H_
#define _UNIFORM01K_H_

#include "RandomFactory.h"

/**
 * @brief A uniform random number generator.
 *
 * This is a uniform generator which draw double between 0 and 1 excluding
 * the end points.
 * This generator is based on an algorithm described by D.E. Knuth, 1981,
 * "Seminumerical Algorithms" 2nd ed., vol.2 of "The Art of Computer
 * Programming" (Reading, MA: Addison-Wesley), §§3.2-3.3.
 * It is addapted from "Numerical Recipes in C".
 */
class Uniform01K : public RandomFactory {
	public: // Constructors and destructor
		/**
		 * @brief Create a Random Number Generator.
		 *
		 * @param seed The seed for the random numbers.
		 */
		Uniform01K(long seed);

		/**
		 * @brief Destroy the generator.
		 */
		~Uniform01K();

	public:
		/**
		 * @brief Set the seed for a new set of random numbers.
		 */
		void setSeed(long seed);

		/**
		 * @brief Get a random number between 0.0 and 1.0 (exclusive of the end point values).
		 */
		double drawNumber() const;

	private:
		static const long MAXNUMBER;
		static const long ZERO;
		static const long MODSEED;
		mutable long _tab[56];
		mutable unsigned int _it1;
		mutable unsigned int _it2;
};
#endif // _UNIFORM01K_H_
