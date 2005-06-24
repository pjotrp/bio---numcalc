//
// File: Parametrizable.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Oct 19 23:06:42 2003
//

#ifndef _PARAMETRIZABLE_H_
#define _PARAMETRIZABLE_H_

// From the STL:
#include <string>
using namespace std;

#include "ParameterList.h"

/**
 * @brief This is the interface for all objects that imply parameters.
 *
 * @see Parameter, ParameterList
 */
class Parametrizable
{
	public:
		Parametrizable() {}
		virtual ~Parametrizable() {}

	public:

		/**
		 * @brief Get all parameters available.
		 *
		 * @return A list with all parameters available.
		 */
		virtual ParameterList getParameters() const = 0;
	
		/**
		 * @brief Get the value for parameter of name 'name'.
		 *
		 * @param name The name of the parameter.
		 * @return the value of parameter <i>name</i>.
		 */
		virtual double getParameterValue(const string & name) const
			throw (ParameterNotFoundException) = 0;

		/**
		 * @brief Set the parameters values to be equals to those of <i>params</i>.
		 *
		 * The list must contain exactly the same parameters (ie same names)
		 * than the parameters available.
		 *
		 * @param parameters A list with all parameters.
		 * @throw ParameterNotFoundException If a some parameter in the list is not in <i>params</i>.
		 * @throw ConstraintException If a value in <i>params</i> does not match the constraint in the
		 * corresponding parameter in the list.
		 */
		virtual void setAllParametersValues(const ParameterList & parameters) 
			throw (ParameterNotFoundException, ConstraintException) = 0;

		/**
		 * @brief Set the value of parameter with name <i>name</i> to be equal to <i>value</i>.
		 *
		 * @param name the name of the parameter to set.
		 * @param value The value of the parameter.
		 * @throw ParameterNotFoundException If no parameter in the list has the name <i>name</i>.
		 * @throw ConstraintException If <i>value</i> does not match the constraint associated to
		 * parameter <i>name</i>.
		 */
		virtual void setParameterValue(const string & name, double value) 
			throw (ParameterNotFoundException, ConstraintException) = 0;

		/**
		 * @brief Update the parameters from <i>params</i>.
		 *
		 * <i>params</i> must be a subset of all parameters available.
		 *
		 * @param parameters A list containing all parameters to update.
		 * @throw ParameterNotFoundException If a some parameter in <i>params</i> is not in the list.
		 * @throw ConstraintException If a value in <i>parameters</i> does not match the constraint in the
		 * corresponding parameter in the list.
		 */
		virtual void setParametersValues(const ParameterList & parameters)
			throw (ParameterNotFoundException, ConstraintException) = 0;

		/**
		 * @brief Update the parameters from <i>parameters</i>.
		 *
		 * Only common parameters with <i>parameters</i> will be updated.
		 *
		 * @param parameters A list of parameters.
		 * @throw ConstraintException If a value in <i>parameters</i> does not match the constraint in the
		 * corresponding parameter in the list.
		 */
		virtual void matchParametersValues(const ParameterList & parameters)
			throw (ConstraintException) = 0;

};

/**
 * @brief A low-level implementation of Parametrizable with void functions.
 *
 * @see Parameter, ParameterList, Parametrizable
 */
class ParametrizableAdapter : public virtual Parametrizable
{
	public:
		ParametrizableAdapter() {}
		virtual ~ParametrizableAdapter() {}

	public:

		/**
		 * @name The Parametrizable interface.
		 *
		 * @{
		 */
		ParameterList getParameters() const { return ParameterList(); }
		double getParameterValue(const string & name) const
			throw (ParameterNotFoundException) { return 0; };
		void setAllParametersValues(const ParameterList & parameters) 
			throw (ParameterNotFoundException, ConstraintException) {}
		void setParameterValue(const string & name, double value) 
			throw (ParameterNotFoundException, ConstraintException) {}
		void setParametersValues(const ParameterList & parameters)
			throw (ParameterNotFoundException, ConstraintException) {}
		void matchParametersValues(const ParameterList & parameters)
			throw (ConstraintException) {};
		/** @} */

};

class AbstractParametrizable : public virtual Parametrizable
{
	protected:
		ParameterList _parameters;
	
	public:
		AbstractParametrizable() {}
		virtual ~AbstractParametrizable() {}

	public:

		ParameterList getParameters() const { return _parameters; }
	
		double getParameterValue(const string & name) const
			throw (ParameterNotFoundException)
		{ return _parameters.getParameter(name) -> getValue(); }

		void setAllParametersValues(const ParameterList & parameters) 
			throw (ParameterNotFoundException, ConstraintException)
		{ _parameters.setAllParametersValues(parameters); }

		void setParameterValue(const string & name, double value) 
			throw (ParameterNotFoundException, ConstraintException)
		{ _parameters.setParameterValue(name, value); }

		void setParametersValues(const ParameterList & parameters)
			throw (ParameterNotFoundException, ConstraintException)
		{ _parameters.setParametersValues(parameters); }

		void matchParametersValues(const ParameterList & parameters)
			throw (ConstraintException)
		{ _parameters.matchParametersValues(parameters); }

};


#endif	//_PARAMETRIZABLE_H_
