//
// File: Functions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Nov  9 23:11:00 2003
//

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "ParameterList.h"
#include "Parametrizable.h"
#include "ParameterExceptions.h"

/**
 * @brief This is the function abstract class.
 */
class Function : public Parametrizable
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
		 * @brief Get the current point at which the function is evaluated.
		 *
		 * @return The set of parameters corresponding to the current point.
		 * @throw Exception If no point is defined.
		 */
		virtual ParameterList getParameters() const throw (Exception) = 0;

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
		virtual double f(const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getValue();
		}
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
		virtual double df(const string & variable, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getFirstOrderDerivative(variable);
		}
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
		 * @param parameters The parameter set to pass to the function.
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




class InfinityFunctionWrapper : public Function
{
	protected:
		Function * _function;
		mutable bool _constraintMatch;
		
	public:
		InfinityFunctionWrapper(Function * function): _function(function), _constraintMatch(false) {}
		virtual ~InfinityFunctionWrapper() {}

	public:

		void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setParameters(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}

		ParameterList getParameters() const throw (Exception) {
			return _function -> getParameters();	
		}

		double getValue() const throw (Exception) {
			return _constraintMatch ? -log(0.) :	_function -> getValue();
		}
		
		double f(const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getValue();
		}
		
		double getParameter(const string & name) const throw (ParameterNotFoundException) {
			return _function -> getParameter(name);
		}
			
		void setAllParametersValues(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setAllParametersValues(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}
		
		void setParameterValue(const string & name, double value) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setParameterValue(name, value);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}
		
		void setParametersValues(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setParametersValues(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}
		
		void matchParametersValues(const ParameterList & parameters) throw (ConstraintException) {
			try {
				_function -> matchParametersValues(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}

};


class InfinityDerivableFirstOrderWrapper : public InfinityFunctionWrapper
{
	public:
		InfinityDerivableFirstOrderWrapper(DerivableFirstOrder * function): InfinityFunctionWrapper(function) {}
		virtual ~InfinityDerivableFirstOrderWrapper() {}

	public:
		
		double getFirstOrderDerivative(const string & variable) const throw (Exception) {
			return _constraintMatch ? -log(0.) :	(dynamic_cast<DerivableFirstOrder *>(_function) -> getFirstOrderDerivative(variable));		
		}
		
		double df(const string & variable, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getFirstOrderDerivative(variable);
		}
};


class InfinityDerivableSecondOrderWrapper : public InfinityDerivableFirstOrderWrapper
{
	public:
		InfinityDerivableSecondOrderWrapper(DerivableFirstOrder * function): InfinityDerivableFirstOrderWrapper(function) {}
		virtual ~InfinityDerivableSecondOrderWrapper() {}

	public:

		double getSecondOrderDerivative(const string & variable) const throw (Exception) {
			return _constraintMatch ? -log(0.) :	(dynamic_cast<DerivableSecondOrder *>(_function) -> getSecondOrderDerivative(variable));					
		}
	
		double d2f(const string & variable, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getSecondOrderDerivative(variable);
		}		

	  double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) {
			return _constraintMatch ? -log(0.) :	(dynamic_cast<DerivableSecondOrder *>(_function) -> getSecondOrderDerivative(variable1, variable2));			
		}
		
		double d2f(const string & variable1, const string & variable2, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getSecondOrderDerivative(variable1, variable2);
		}

};

#endif	//_FUNCTIONS_H_
