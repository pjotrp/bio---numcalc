//
// File: ParameterList.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Oct 15 18:17:29 2003
//

#ifndef _PARAMETERLIST_H_
#define _PARAMETERLIST_H_

#include "Parameter.h"

// From STL:
#include <vector>
#include <string>
#include <iostream>

using namespace std;

/**
 * @brief The parameter list object.
 * This is a vector of Parameter with a few additional methods, mainly for giving
 * name access.
 */
class ParameterList : public vector<Parameter *>
{
	public: // Class constructor and destructor:
		
		/**
		 * @brief Build a new ParameterList object.
		 */
		ParameterList();
	
		ParameterList(const ParameterList & pl);
		
		ParameterList & operator=(const ParameterList & pl);
	
		virtual ~ParameterList();
	
	public:
		
		/**
		 * @brief Get the parameter with name <i>name</i>.
		 *
		 * @param name The name of the parameter to look for.
		 * @return A const pointer toward the parameter with name <i>name</i>, or NULL if not found.
		 */
		virtual const Parameter * getParameter(const string & name) const;
	
		/**
		 * @brief Get the parameter with name <i>name</i>.
		 *
		 * @param name The name of the parameter to look for.
		 * @return A pointer toward the parameter with name <i>name</i>, or NULL if not found.
		 */
		virtual Parameter * getParameter(const string & name);
	
		/**
		 * @brief Get all parameter names in the list.
		 *
		 * @return A vector with all names in the same order as the parameters in the list.
		 */
		virtual vector<string> getParameterNames() const;
	
		/**
		 * @brief Add a new parameter at the end of the list.
		 *
		 * @param param The parameter to add to the list.
		 */
		virtual void addParameter(const Parameter & param) throw (ParameterException);
		
		/**
		 * @brief Add new parameters at the end of the list.
		 *
		 * @param params The parameter list containing the new paramters to add to the list.
		 */
		virtual void addParameters(const ParameterList & params) throw (ParameterException);
		
		/**
		 * @brief Set the value of parameter with name <i>name</i> to be equal to <i>value</i>.
		 *
		 * @param name the name of the parameter to set.
		 * @param value The value of the parameter.
		 */
		virtual void setParameterValue(const string & name, double value) 
			throw (ParameterNotFoundException, ConstraintException);

		/**
		 * @brief Set the parameters to be equals to <i>params</i>.
		 *
		 * The list must contain exactly the same parameters (ie same names)
		 * than the parameters available.
		 *
		 * @param params A list with all parameters.
		 * @see setParameters(), matchParameters();
		 */
		virtual void setAllParametersValues(const ParameterList & params)
			throw (ParameterNotFoundException, ConstraintException);

		/**
		 * @brief Update the parameters from <i>params</i>.
		 *
		 * <i>params</i> must be a subset of all parameters available.
		 *
		 * @param params A list containing all parameters to update.
		 * @see setAllParameters(), matchParameters()
		 */
		virtual void setParametersValues(const ParameterList & params)
			throw (ParameterNotFoundException, ConstraintException);

		/**
		 * @brief Update the parameters from <i>params</i>.
		 *
		 * Only common parameters with <i>params</i> will be updated.
		 *
		 * @param params A list of parameters.
		 * @see setParameters(), setAllParameters()
		 */
		virtual void matchParametersValues(const ParameterList & params)
			throw (ConstraintException);

		/**
		 * @brief Set the parameters to be equals to <i>params</i>.
		 *
		 * The list must contain exactly the same parameters (ie same names)
		 * than the parameters available.
		 *
		 * @param params A list with all parameters.
		 * @see setParameters(), matchParameters();
		 */
		virtual void setAllParameters(const ParameterList & params)
			throw (ParameterNotFoundException);

		/**
		 * @brief Update the parameters from <i>params</i>.
		 *
		 * <i>params</i> must be a subset of all parameters available.
		 *
		 * @param params A list containing all parameters to update.
		 * @see setAllParameters(), matchParameters()
		 */
		virtual void setParameters(const ParameterList & params)
			throw (ParameterNotFoundException);

		/**
		 * @brief Update the parameters from <i>params</i>.
		 *
		 * Only common parameters with <i>params</i> will be updated.
		 *
		 * @param params A list of parameters.
		 * @see setParameters(), setAllParameters()
		 */
		virtual void matchParameters(const ParameterList & params);

		/**
		 * @brief Delete a parameter from the list.
		 *
		 * @param name The name of the parameter to delete from the list.
		 */
		virtual void deleteParameter(const string & name) throw (ParameterNotFoundException);
		
		/**
		 * @brief Print all parameters.
		 */
		virtual void printParameters(ostream & out) const;
		
		/**
		 * @brief Reset the list: delete all parameters.
		 */
		virtual void reset();
};


#endif	//_PARAMETERLIST_H_
