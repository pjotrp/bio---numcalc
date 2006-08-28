//
// File: AbstractNumericalDerivative.h
// Created by: Julien Dutheil
// Created on: Thu Aug 17 15:00 2006
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _ABSTRACTNUMERICALDERIVATIVE_H_
#define _ABSTRACTNUMERICALDERIVATIVE_H_

#include "Functions.h"

//From the STL:
#include <map>
#include <vector>
#include <string>

/**
 * @brief Numerical derivative function wrapper, partial implementation.
 *
 * This class provides a wrapper for Function object, implementing the DerivableSecondOrder interface
 * (Although no cross derivative is implemented for now).
 * Derivations of this class can be used as full DerivableSecondOrder objects, with derivative functions.
 * 
 * Three kinds of constructors are provided: one with a Function object, another with a DerivableFirstOrder object, and one with a DerivableSecondOrder object.
 * In the first case, all derivatives will be computed analytically.
 * In the second case, first order derivative will be computed analytically only if no appropriate analytical derivative is available, second order derivative will always be computed numerically.
 * In the last case, first and second order derivative will be computed analytically only if no appropriate analytical derivative is available.
 */
class AbstractNumericalDerivative: public DerivableSecondOrder
{
  protected:
    Function * _function;
    DerivableFirstOrder * _function1;
    DerivableSecondOrder * _function2;
    double _h;
    vector<string> _variables;
    mutable map<string,double> _der1;
    mutable map<string,double> _der2;
    
	public:
		AbstractNumericalDerivative (Function * function):
      _function(function), _function1(NULL), _function2(NULL), _h(0.0001) {}
		AbstractNumericalDerivative (DerivableFirstOrder * function):
      _function(function), _function1(function), _function2(NULL), _h(0.0001) {}
	  AbstractNumericalDerivative (DerivableSecondOrder * function):
      _function(function), _function1(function), _function2(function), _h(0.0001) {}
		virtual ~AbstractNumericalDerivative() {}

  public:
    /**
     * @brief Set the interval value used in numerical approximation.
     *
     * Default value is 0.0001.
     *
     * @param h Interval value.
     */
    void setInterval(double h) { _h = h; }
    
    /**
     * @brief Set the list of parameters to derivate.
     *
     * @param variables A list of all parameter names.
     */
    void setParametersToDerivate(const vector<string> & variables)
    {
      _variables = variables;
      _der1.clear();
      _der2.clear();
    }
    
    /**
     * @name The DerivableFirstOrder interface
     *
     * @{
     */
    double getFirstOrderDerivative(const string & variable) const
      throw (Exception)
    {
      if(_function1 != NULL)
      {
        try {
          return _function1->getFirstOrderDerivative(variable);
        } catch(Exception & e) {}
      }    
      map<string,double>::iterator it = _der1.find(variable);
      if(it != _der1.end()) return it->second;
      else throw Exception("First order derivative not computed for variable " + variable + "."); 
    }
    /** @} */
		
    /**
     * @name The DerivableSecondOrder interface
     *
     * @{
     */
    double getSecondOrderDerivative(const string & variable) const
      throw (Exception)
    {
      if(_function2 != NULL)
      {
        try {
          return _function2->getSecondOrderDerivative(variable);
        } catch(Exception & e) {}
      }    
      map<string,double>::iterator it = _der2.find(variable);
      if(it != _der2.end()) return it->second;
      else throw Exception("Second order derivative not computed for variable " + variable + "."); 
    }

		double getSecondOrderDerivative(const string & variable1, const string & variable2) const
      throw (Exception)
    {
      throw Exception("Unimplemented cross derivative.");
    }
	  
		ParameterList getParameters() const throw (Exception)
    {
			return _function->getParameters();	
		}
    /** @} */
		
    
    /**
     * @name The Function interface
     *
     * @{
     */
		double getParameterValue(const string & name) const
      throw (ParameterNotFoundException)
    {
			return _function->getParameterValue(name);
		}
			
		void setAllParametersValues(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
			_function->setAllParametersValues(parameters);
		}
		
		void setParameterValue(const string & name, double value)
      throw (ParameterNotFoundException, ConstraintException)
    {
			_function->setParameterValue(name, value);
		}
		
		void setParametersValues(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
			_function->setParametersValues(parameters);
		}
		
		void matchParametersValues(const ParameterList & parameters)
      throw (ConstraintException)
    {
			_function -> matchParametersValues(parameters);
		}
    /** @} */
    
};

#endif //_ABSTRACTNUMERICALDERIVATIVE_H_

