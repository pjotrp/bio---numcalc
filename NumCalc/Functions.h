//
// File: Functions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Nov  9 23:11:00 2003
//

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "ParameterList.h"

/**
 * @brief This is the function interface.
 */
class Function
{
	public:
		Function() {}
		virtual ~Function() {}

	public:
		/**
		 * @brief Get the value of the function according to a given set of parameters.
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 */
		virtual double f(const ParameterList & parameters) const = 0;
};

/**
 * @brief This is the interface for first order derivable functions.
 */
class DerivableFirstOrder : public Function
{
	public:
		DerivableFirstOrder() {}
		virtual ~DerivableFirstOrder() {}

	public:

		/**
		 * @brief Get the value of the first derivative of the function
		 * according to a given set of parameters.
		 *
		 * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{df}{dx} @f$.
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 */
		virtual double df(const string & variable, const ParameterList & parameters) const = 0;		
};

/**
 * @brief This is the interface for second order derivable functions.
 */
class DerivableSecondOrder : public DerivableFirstOrder
{
	public:
		DerivableSecondOrder() {}
		virtual ~DerivableSecondOrder() {}

	public:

		/**
		 * @brief Get the value of the second derivative of the function
		 * according to a given set of parameters.
		 *
		 * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x^2} @f$.
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 */
		virtual double d2f(const string & variable, const ParameterList & parameters) const = 0;		

		/**
		 * @brief Get the value of the cross derivative of the function
		 * according to a given set of parameters.
		 *
		 * @param variable1  The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
		 * @param variable2  The name of the @f$ y @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 */
		virtual double d2f(const string & variable1, const string & variable2, const ParameterList & parameters) const = 0;		
};

#endif	//_FUNCTIONS_H_
