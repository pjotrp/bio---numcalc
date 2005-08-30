//
// File: Parameter.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Oct 15 15:40:47 2003
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
