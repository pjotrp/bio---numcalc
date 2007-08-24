//
// File: Parametrizable.h
// Created by: Julien Dutheil
// Created on: Sun Oct 19 23:06:42 2003
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

#ifndef _PARAMETRIZABLE_H_
#define _PARAMETRIZABLE_H_

// From Utils:
#include <Utils/Clonable.h>

// From the STL:
#include <string>
using namespace std;

#include "ParameterList.h"

/**
 * @brief This is the interface for all objects that imply parameters.
 *
 * @see Parameter, ParameterList
 */
class Parametrizable:
  public virtual Clonable
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
     * @brief Get the parameter with specified name.
     *
     * @param name The name of the parameter to look for.
     * @return The parameter with given name.
     * @throw ParameterNotFoundException if no parameter with this name is found.
     */
    virtual Parameter getParameter(const string & name) const throw (ParameterNotFoundException) = 0;
	
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
 * @brief A low-level implementation of the Parametrizable interface with void functions.
 *
 * @see Parameter, ParameterList, Parametrizable
 */
class ParametrizableAdapter:
  public virtual Parametrizable
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
    Parameter getParameter(const string & name) const throw (ParameterNotFoundException) { return Parameter(); }
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

/**
 * @brief A partial implementation of the Parametrizable interface.
 *
 * Parameters are stored in a protected ParameterList object.
 *
 * The abstract fireParameterChanged() method is provided so that the derived class
 * know when a parameter has changed, and can be updated.
 * All methods call the corresponding method in ParameterList and then call the
 * fireParameterChanged() method.
 */
class AbstractParametrizable:
  public virtual Parametrizable
{
	protected:
		mutable ParameterList _parameters;
	
	public:
		AbstractParametrizable() {}

		virtual ~AbstractParametrizable() {}

	public:

		ParameterList getParameters() const { return _parameters; }
    
    Parameter getParameter(const string & name) const throw (ParameterNotFoundException)
    {
      const Parameter * p = _parameters.getParameter(name);
      if(p) return *p;
      else throw ParameterNotFoundException("AbstractParametrizable::getParameter.", name);
    }
	
		double getParameterValue(const string & name) const
			throw (ParameterNotFoundException)
		{ 
			return _parameters.getParameter(name)->getValue();
		}

		void setAllParametersValues(const ParameterList & parameters) 
			throw (ParameterNotFoundException, ConstraintException)
		{
			_parameters.setAllParametersValues(parameters);
			fireParameterChanged(parameters);
		}

		void setParameterValue(const string & name, double value) 
			throw (ParameterNotFoundException, ConstraintException)
		{ 
			_parameters.setParameterValue(name, value);
			fireParameterChanged(_parameters.subList(name));
		}

		void setParametersValues(const ParameterList & parameters)
			throw (ParameterNotFoundException, ConstraintException)
		{ 
			_parameters.setParametersValues(parameters);
			fireParameterChanged(parameters);
		}

		void matchParametersValues(const ParameterList & parameters)
			throw (ConstraintException)
		{ 
			_parameters.matchParametersValues(parameters);
			fireParameterChanged(parameters);
		}

		/**
		 * @brief Notify the class when one or several parameters have changed.
		 *
		 * @param parameters A ParameterList object with parameters that changed.
		 */
		virtual void fireParameterChanged(const ParameterList & parameters) = 0;

};

#endif	//_PARAMETRIZABLE_H_

