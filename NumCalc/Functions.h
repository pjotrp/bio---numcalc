//
// File: Functions.h
// Created by: Julien Dutheil
// Created on: Sun Nov  9 23:11:00 2003
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

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "ParameterList.h"
#include "Parametrizable.h"
#include "ParameterExceptions.h"

// From Utils:
#include <Utils/Clonable.h>

// From the STL:
#include <cmath>
using namespace std;

/**
 * @brief This is the function abstract class.
 *
 * This class provides the interface for function objet
 * and a default implementation of the f() function.
 *
 * The f() function sends the value of the function according to a
 * given set of parameters.
 *
 * However for complexe function like likelihood for instance,
 * computing the function value takes some time, and one do not want
 * to perform several times the computation for an identical set of 
 * parameters.
 * The setParameters() method hence allows to set the parameter value
 * for which the function is to be computed, perform the computation
 * and store the results.
 * The getValue() methods send the result of the computation.
 * One may hence access to the result of the computation by calling the
 * getvalue() method without re-computating the function.
 * The f(parameters) function is a shortcut for
 * @code
 * setParameters(parameters);
 * return getValue();
 * @endcode
 * for convinience.
 *
 * @see Parameter, ParameterList
 */
class Function :
  public virtual Parametrizable,
  public virtual Clonable
{    
  public:
    Function() {}
    virtual ~Function() {}

  public:

    /**
     * @brief Set the point where the function must be computed.
     *
     * @param parameters The parameter set to pass to the function.
     */
    virtual void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) = 0;

    /**
     * @brief Get the value of the function at the current point.
     *
     * @return The value of the function.
     * @throw Exception If no point is specified or if an error occured.
     */
    virtual double getValue() const throw (Exception) = 0;
    
    /**
     * @brief Get the value of the function according to a given set of parameters.
     * 
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double f(const ParameterList & parameters) throw (Exception)
    {
      setParameters(parameters);
      return getValue();
    }
};

/**
 * @brief This is the abstract class for first order derivable functions.
 *
 * This class adds the getFirstOrderDerivative() and df() shortcut functions.
 */
class DerivableFirstOrder : public virtual Function
{
  public:
    DerivableFirstOrder() {}
    virtual ~DerivableFirstOrder() {}

  public:
    
    /**
     * @brief Get the derivative of the function at the current point.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{df}{dx} @f$.
     * @return The value of the function.
     * @throw Exception If no point is specified or if an error occured.
     */
    virtual double getFirstOrderDerivative(const string & variable) const throw (Exception) = 0;
    
    /**
     * @brief Get the value of the first derivative of the function
     * according to a given set of parameters.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{df}{dx} @f$.
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double df(const string & variable, const ParameterList & parameters) throw (Exception)
    {
      setParameters(parameters);
      return getFirstOrderDerivative(variable);
    }
};

/**
 * @brief This is the abstract class for second order derivable functions.
 * 
 * This class adds the getSecondOrderDerivative() and d2f() shortcut functions.
 * Cross derivative functions are also provided.
 */
class DerivableSecondOrder : public virtual DerivableFirstOrder
{
  public:
    DerivableSecondOrder() {}
    virtual ~DerivableSecondOrder() {}

  public:

    /**
     * @brief Get the second order derivative of the function at the current point.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x^2} @f$.
     * @return The value of the function.
     * @throw Exception If no point is specified or if an error occured.
     */
    virtual double getSecondOrderDerivative(const string & variable) const throw (Exception) = 0;
  
    /**
     * @brief Get the value of the second order derivative of the function
     * according to a given set of parameters.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x^2} @f$.
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double d2f(const string & variable, const ParameterList & parameters) throw (Exception) {
      setParameters(parameters);
      return getSecondOrderDerivative(variable);
    }    

    /**
     * @brief Get the value of the cross derivative of the function
     * according to a given set of parameters.
     *
     * @param variable1  The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @param variable2  The name of the @f$ y @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) = 0;  
    
    /**
     * @brief Get the value of the cross derivative of the function
     * according to a given set of parameters.
     *
     * @param variable1  The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @param variable2  The name of the @f$ y @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double d2f(const string & variable1, const string & variable2, const ParameterList & parameters) throw (Exception) {
      setParameters(parameters);
      return getSecondOrderDerivative(variable1, variable2);
    }
};



/**
 * @brief Wrapper class for optimization under constraints.
 *
 * Catch any ConstraintException thrown and send +inf.
 */
class InfinityFunctionWrapper : public virtual Function
{
  protected:
    Function * _function;
    mutable bool _constraintMatch;
    
  public:
    InfinityFunctionWrapper(Function * function):
      _function(function),
      _constraintMatch(false) {}
    virtual ~InfinityFunctionWrapper() {}

#if defined(VIRTUAL_COV)
    InfinityFunctionWrapper * clone() const { return new InfinityFunctionWrapper(*this); }
#else
    Clonable * clone() const { return new InfinityFunctionWrapper(*this); }
#endif

  public:

    void setParameters(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      try {
        _function->setParameters(parameters);
        _constraintMatch = false;
      } catch(ConstraintException ce) {
        _constraintMatch = true;
      }
    }

    ParameterList getParameters() const throw (Exception)
    {
      return _function->getParameters();  
    }

    double getValue() const throw (Exception)
    {
      return _constraintMatch ? -log(0.) :  _function->getValue();
    }
    
    double f(const ParameterList & parameters) throw (Exception)
    {
      setParameters(parameters);
      return getValue();
    }
    
    double getParameterValue(const string & name) const throw (ParameterNotFoundException)
    {
      return _function->getParameterValue(name);
    }
      
    void setAllParametersValues(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      try {
        _function->setAllParametersValues(parameters);
        _constraintMatch = false;
      } catch(ConstraintException ce) {
        _constraintMatch = true;
      }
    }
    
    void setParameterValue(const string & name, double value)
      throw (ParameterNotFoundException, ConstraintException)
    {
      try {
        _function->setParameterValue(name, value);
        _constraintMatch = false;
      } catch(ConstraintException ce) {
        _constraintMatch = true;
      }
    }
    
    void setParametersValues(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      try {
        _function->setParametersValues(parameters);
        _constraintMatch = false;
      } catch(ConstraintException ce) {
        _constraintMatch = true;
      }
    }
    
    void matchParametersValues(const ParameterList & parameters)
      throw (ConstraintException)
    {
      try {
        _function->matchParametersValues(parameters);
        _constraintMatch = false;
      } catch(ConstraintException ce) {
        _constraintMatch = true;
      }
    }

};

/**
 * @brief Wrapper class for optimization under constraints.
 *
 * Catch any ConstraintException thrown and send +inf.
 */
class InfinityDerivableFirstOrderWrapper : public virtual InfinityFunctionWrapper
{
  public:
    InfinityDerivableFirstOrderWrapper(DerivableFirstOrder * function): InfinityFunctionWrapper(function) {}
    virtual ~InfinityDerivableFirstOrderWrapper() {}
    
#if defined(VIRTUAL_COV)
    InfinityDerivableFirstOrderWrapper * clone() const { return new InfinityDerivableFirstOrderWrapper(*this); }
#else
    Clonable * clone() const { return new InfinityDerivableFirstOrderWrapper(*this); }
#endif

  public:
    
    double getFirstOrderDerivative(const string & variable) const throw (Exception) {
      return _constraintMatch ? -log(0.) :  (dynamic_cast<DerivableFirstOrder *>(_function)->getFirstOrderDerivative(variable));    
    }
    
    double df(const string & variable, const ParameterList & parameters) throw (Exception) {
      setParameters(parameters);
      return getFirstOrderDerivative(variable);
    }
};

/**
 * @brief Wrapper class for optimization under constraints.
 *
 * Catch any ConstraintException thrown and send +inf.
 */
class InfinityDerivableSecondOrderWrapper : public virtual InfinityDerivableFirstOrderWrapper
{
  public:
    InfinityDerivableSecondOrderWrapper(DerivableFirstOrder * function):
      InfinityFunctionWrapper(function),
      InfinityDerivableFirstOrderWrapper(function) {}
    virtual ~InfinityDerivableSecondOrderWrapper() {}

#if defined(VIRTUAL_COV)
    InfinityDerivableSecondOrderWrapper * clone() const { return new InfinityDerivableSecondOrderWrapper(*this); }
#else
    Clonable * clone() const { return new InfinityDerivableSecondOrderWrapper(*this); }
#endif

  public:

    double getSecondOrderDerivative(const string & variable) const throw (Exception) {
      return _constraintMatch ? -log(0.) :  (dynamic_cast<DerivableSecondOrder *>(_function)->getSecondOrderDerivative(variable));          
    }
  
    double d2f(const string & variable, const ParameterList & parameters) throw (Exception) {
      setParameters(parameters);
      return getSecondOrderDerivative(variable);
    }    

    double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) {
      return _constraintMatch ? -log(0.) :  (dynamic_cast<DerivableSecondOrder *>(_function) -> getSecondOrderDerivative(variable1, variable2));      
    }
    
    double d2f(const string & variable1, const string & variable2, const ParameterList & parameters) throw (Exception) {
      setParameters(parameters);
      return getSecondOrderDerivative(variable1, variable2);
    }

};

#endif  //_FUNCTIONS_H_

