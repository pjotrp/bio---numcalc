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
#include "Matrix.h"

//From the STL:
#include <map>
#include <vector>
#include <string>

using namespace std;

namespace bpp
{

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
class AbstractNumericalDerivative:
  public DerivableSecondOrder,
  public FunctionWrapper
{
  protected:
    DerivableFirstOrder *_function1;
    DerivableSecondOrder *_function2;
    double _h;
    vector<string> _variables;
    mutable map<string, unsigned int> _index; //Store positions in array corresponding to variable names.
    vector<double> _der1;
    vector<double> _der2;
    RowMatrix<double> _crossDer2;
    bool _computeD1, _computeD2, _computeCrossD2;
    
	public:
		AbstractNumericalDerivative (Function * function):
      FunctionWrapper(function), _function1(NULL), _function2(NULL), _h(0.0001), _computeD1(true), _computeD2(true), _computeCrossD2(false) {}
		AbstractNumericalDerivative (DerivableFirstOrder * function):
      FunctionWrapper(function), _function1(function), _function2(NULL), _h(0.0001), _computeD1(true), _computeD2(true), _computeCrossD2(false) {}
	  AbstractNumericalDerivative (DerivableSecondOrder * function):
      FunctionWrapper(function), _function1(function), _function2(function), _h(0.0001), _computeD1(true), _computeD2(true), _computeCrossD2(false) {}
		virtual ~AbstractNumericalDerivative() {}

#ifndef NO_VIRTUAL_COV
    AbstractNumericalDerivative*
#else
    Clonable*
#endif
    clone() const = 0;

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
      _index.clear();
      for(unsigned int i = 0; i < _variables.size(); i++)
        _index[_variables[i]] = i;
      _der1.resize(_variables.size());
      _der2.resize(_variables.size());
      _crossDer2.resize(_variables.size(), _variables.size());
    }
    
    /**
     * @name The DerivableFirstOrder interface
     *
     * @{
     */
    void enableFirstOrderDerivatives(bool yn) { _computeD1 = yn; }
    bool enableFirstOrderDerivatives() const { return _computeD1; }
    
    double getFirstOrderDerivative(const string & variable) const
      throw (Exception)
    {
      if(_function1 != NULL)
      {
        try
        {
          return _function1->getFirstOrderDerivative(variable);
        }
        catch(Exception & e) {}
      }    
      map<string, unsigned int>::iterator it = _index.find(variable);
      if(_computeD1 && it != _index.end()) return _der1[it->second];
      else throw Exception("First order derivative not computed for variable " + variable + "."); 
    }
    /** @} */
		
    /**
     * @name The DerivableSecondOrder interface
     *
     * @{
     */

    void enableSecondOrderDerivatives(bool yn) { _computeD2 = yn; }
    bool enableSecondOrderDerivatives() const { return _computeD2; }

    double getSecondOrderDerivative(const string & variable) const
      throw (Exception)
    {
      if(_function2 != NULL)
      {
        try
        {
          return _function2->getSecondOrderDerivative(variable);
        }
        catch(Exception & e) {}
      }    
      map<string, unsigned int>::iterator it = _index.find(variable);
      if(_computeD2 && it != _index.end()) return _der2[it->second];
      else throw Exception("Second order derivative not computed for variable " + variable + "."); 
    }

		double getSecondOrderDerivative(const string & variable1, const string & variable2) const
      throw (Exception)
    {
      //throw Exception("Unimplemented cross derivative.");
      if(_function2 != NULL)
      {
        try
        {
          return _function2->getSecondOrderDerivative(variable1, variable2);
        }
        catch(Exception & e) {}
      }    
      map<string, unsigned int>::iterator it1 = _index.find(variable1);
      map<string, unsigned int>::iterator it2 = _index.find(variable2);
      if(_computeCrossD2 && it1 != _index.end() && it2 != _index.end()) return _crossDer2(it1->second, it2->second);
      else throw Exception("Cross second order derivative not computed for variables " + variable1 + " and " + variable2 + "."); 
    }
    /** @} */
	   
    /**
     * @name The Parametrizable interface.
     *
     * @{
     */
    double f(const ParameterList & parameters) throw (Exception)
    {
      setParameters(parameters);
      return getValue();
    }
    void setParameters(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      _function->setParameters(parameters);
      updateDerivatives(parameters);
    }
    void setAllParametersValues(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      _function->setAllParametersValues(parameters);
      updateDerivatives(parameters);
    }
    
    void setParameterValue(const string & name, double value)
      throw (ParameterNotFoundException, ConstraintException)
    {
      _function->setParameterValue(name, value);
      updateDerivatives(_function->getParameters().subList(name));
    }
    
    void setParametersValues(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      _function->setParametersValues(parameters);
      updateDerivatives(parameters);
    }
    
    void matchParametersValues(const ParameterList & parameters)
      throw (ConstraintException)
    {
      _function->matchParametersValues(parameters);
      updateDerivatives(parameters);
    }
    /** @} */

    void enableSecondOrderCrossDerivatives(bool yn) { _computeCrossD2 = yn; }
    bool enableSecondOrderCrossDerivatives() const { return _computeCrossD2; }

  protected:
    virtual void updateDerivatives(const ParameterList & parameters) = 0;
    
};

} //end of namespace bpp.

#endif //_ABSTRACTNUMERICALDERIVATIVE_H_

