//
// File: Parameter.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Oct 15 15:40:47 2003
//

#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include "ParameterExceptions.h"
#include "Constraints.h"

// From the STL:
#include <string>

//From Utils:
#include <Utils/Clonable.h>

/**
 * @brief This class is designed to facilitate the manipulation of parameters.
 *
 * A parameter contains a <i>value</i> stored as a double.
 * It also contains a <i>name</i> and optionaly a Constraint.
 * The constraint object allows to apply restriction on the value of the parameter,
 * for instance positive number, or a particular interval and so on.
 *
 * @see ParameterList, Parametrizable
 */
class Parameter: public Clonable
{
	protected:
		string name;					//Parameter name
		double value;					//Parameter value
		const Constraint * constraint; 	//A constraint on the value
	
	public: // Class constructors and destructors:
		
		/**
		 * @brief Build a new parameter.
		 *
		 * @param name       The parameter name.
		 * @param value      The parameter value.
		 * @param constraint An optional pointer toward a constraint Object.
		 * @throw ConstraintException If the parameter value does not match the contraint.
		 */
		Parameter(const string & name = "", double value = 0, const Constraint * constraint = NULL)
		throw (ConstraintException);

		//Copy constructor:
		/**
		 * @brief Copy constructor.
		 */
		Parameter(const Parameter & param);
		/**
		 * @brief Assignment operator.
		 */
		virtual Parameter & operator=(const Parameter & param);
	
		//Destructor:
		virtual ~Parameter();
		
	public: // The clonable interface is implemented here:
		
		Clonable * clone() const;

	public:
		/**
		 * @brief Set the name of this parameter.
		 *
		 * @param name the new parameter name.
		 */
		virtual void setName(const string & name);
	
		/**
		 * @brief Set the value of this parameter.
		 *
		 * @param value the new parameter value.
		 */
		virtual void setValue(double value) throw (ConstraintException);
	
		/**
		 * @brief Get the name of this parameter.
		 *
		 * @return The parameter name.
		 */
		virtual string getName() const;
	
		/**
		 * @brief Get the value of this parameter.
		 *
		 * @return The parameter value.
		 */
		virtual double getValue() const;
		
		/**
		 * @brief Return the constraint associated to this parameter if there is one.
		 *
		 * @return A pointer toward the constraint, or NULL if there is no constraint.
		 */
		virtual const Constraint * getConstraint() const;

		/**
		 * @brief Tells if this parameter has a constraint.
		 *
		 * @return True if this parameter has a contraint.
		 */
		virtual bool hasConstraint() const;
		
		/**
		 * @brief Remove the constraint assoviated to this parameter.
		 *
		 * @return A pointer toward the formerly used contraint.
		 */
		virtual const Constraint * removeConstraint();
	
	public:
		static IncludingPositiveReal R_PLUS;
		static ExcludingPositiveReal R_PLUS_STAR;
};

#endif	//_PARAMETER_H_
