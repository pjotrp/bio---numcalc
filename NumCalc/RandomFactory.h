/*
 * File RandomFactory.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
 */

// Secured inclusion of header's file
#ifndef _RANDOMFACTORY_H_
#define _RANDOMFACTORY_H_

/**
 * @brief This is the interface for the Random Number Generators.
 *
 * A Random Number Generator draw numbers between two end points. Each
 * number is taken from a given statistic distribution.
 */
class RandomFactory {
	public:
		virtual ~RandomFactory() {}
	public:
		/**
		 * @brief Set the generator's seed.
		 */
		virtual void setSeed(long seed) = 0;

		/**
		 * @brief Return a random number.
		 */
		virtual double drawNumber() const = 0;
};

#endif // _RANDOMFACTORY_H_
