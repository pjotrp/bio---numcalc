//
// File: ParameterList.h
// Created by: Julien Dutheil
// Created on: Wed Oct 15 18:17:29 2003
//

/*
Copyright or © or Copr. CNRS, (November 19, 2004)

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
 * 
 * This is a vector of Parameter with a few additional methods, mainly for giving
 * name access.
 */
class ParameterList : public virtual vector<Parameter *>
{
	public:
		
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
		 * @brief Get given parameters as a sublist.
		 *
		 * @param names Name of the parameters to be included in the list.
		 * @return A list with all parameters specified.
		 */
		virtual ParameterList subList(const vector<string> & names) const ;
		
		/**
		 * @brief Get given parameter as a sublist.
		 *
		 * @param name Name of the parameter to be included in the list.
		 * @return A list with the parameter specified.
		 */
		virtual ParameterList subList(const string & name) const;

		/**
		 * @brief Get given parameters as a sublist.
		 *
		 * @param parameters Positions of the parameters to be included in the list.
		 * @return A list with all parameters specified.
		 */
		virtual ParameterList subList(vector<unsigned int> parameters) const;

		/**
		 * @brief Get given parameter as a sublist.
		 *
		 * @param parameter Position of the parameters to be included in the list.
		 * @return A list with the parameter specified.
		 */
		virtual ParameterList subList(unsigned int parameter) const;

		/**
		 * @brief Get the sublist containing all common parameter between this list and pl.
		 *
		 * @param params The list to compare to.
		 * @return A list with all common parameters.
		 */
		virtual ParameterList getCommonParametersWith(const ParameterList & params) const;
	
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

