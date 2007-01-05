//
// File: AbstractOptimizer.h
// Created by: Julien Dutheil
// Created on: Mon Dec 22 12:18:09 2003
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

#ifndef _ABSTRACTOPTIMIZER_H_
#define _ABSTRACTOPTIMIZER_H_

#include "Optimizer.h"

/**
 * @brief Partial implementation of the Optimizer interface.
 */
class AbstractOptimizer: public Optimizer
{
	protected:
		
		/**
		 * @brief The function to optimize.
		 */
		Function * _function;
	
		/**
		 * @brief The parameters that will be optimized.
		 */
		ParameterList _parameters;
	
		/**
		 * @brief The message handler.
		 */
		ostream * _messageHandler;
	
		/**
		 * @brief The profiler.
		 */
		ostream * _profiler;
		
		/**
		 * @brief The constraint policy.
		 *
		 * Must be one the following:
		 * - CONSTRAINTS_KEEP: keep the constraint associated to the parameters (default).
		 * - CONSTRAINTS_IGNORE: remove all constraints.
		 * - CONSTRAINTS_AUTO: use AutoParameters to deal with constraints.
		 *
		 * @see AutoParameter
		 */ 		 
		string _constraintPolicy;

		/**
		 * @brief Tell if the tolerance level has been reached.
		 *
		 * This field is initilaised by the init() method, maintained by the
		 * step() method and used in the optimize() method.
		 */
		bool _tolIsReached;
		
		/**
		 * @brief The stoping condition to use while optimizing.
		 */
		OptimizationStopCondition * _stopCondition;
		
		/**
		 * @brief The default stoping condition to use while optimizing.
		 */
		OptimizationStopCondition * _defaultStopCondition;

		/**
		 * @brief The maximum number of function evaluations allowed.
		 */
		int _nbEvalMax;
		
		/**
		 * @brief The current number of function evaluations achieved.
		 */
		int _nbEval;

		/**
		 * @brief State of the verbose mode: > 0 = enabled.
		 *
		 * This may not be used by the Optimizer.
		 */
		unsigned int _verbose;


	public:
		AbstractOptimizer(Function * function = NULL);

    AbstractOptimizer(const AbstractOptimizer & opt);
    
    AbstractOptimizer & operator=(const AbstractOptimizer & opt);

		virtual ~AbstractOptimizer()
    {
      delete _stopCondition;
      delete _defaultStopCondition;
    }
	
	public:
		
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */
		void init(const ParameterList & params) throw (Exception);
		ParameterList getParameters() const { return _parameters; }
		void setFunction(Function * function)
    { 
      _function = function;
      if(function != NULL) _stopCondition->init();
    }
		const Function * getFunction() const { return _function; }
		Function * getFunction() { return _function; }
		double getFunctionValue() const throw (NullPointerException)
		{
			if(_function == NULL) throw NullPointerException("AbstractOptimizer::getFunctionValue. No function associated to this optimizer.");
			return _function -> getValue();
		}
		void setMessageHandler(ostream * mh) { _messageHandler = mh; }
		void setProfiler(ostream * profiler) { _profiler = profiler; }
		int getNumberOfEvaluations() const { return _nbEval; }
		void setStopCondition(const OptimizationStopCondition & stopCondition)
    {
      _stopCondition = stopCondition.clone();
    }
		OptimizationStopCondition * getStopCondition() { return _stopCondition; }
		const OptimizationStopCondition * getStopCondition() const { return _stopCondition; }
		OptimizationStopCondition * getDefaultStopCondition() { return _defaultStopCondition; }
		const OptimizationStopCondition * getDefaultStopCondition() const { return _defaultStopCondition; }
		bool isToleranceReached() const { return _tolIsReached; }
		bool isMaximumNumberOfEvaluationsReached() const { return _nbEvalMax >= _nbEvalMax; }
		void setMaximumNumberOfEvaluations(unsigned int max) { _nbEvalMax = max; }
		void setVerbose(unsigned int v) { _verbose = v; }
		unsigned int getVerbose() const { return _verbose; }
		void setConstraintPolicy(const string & constraintPolicy) { _constraintPolicy = constraintPolicy; }
		string getConstraintPolicy() const { return _constraintPolicy; }
		/** @} */
	
	protected:
		
		/**
		 * @name Inner utilitary functions
		 *
		 * @{
		 */
		
		/**
		 * @brief Build a list of AutoParameter instead of Parameter.
		 */
		void autoParameter();
	
		/**
		 * @brief Remove the constraints of all the arguments.
		 */
		void ignoreConstraints();
	
		/**
		 * @brief Print to the profile if there is one.
		 *
		 * @param v The double value to print.
		 */
		void profile(double v);
	
		/**
		 * @brief Print to the profile if there is one.
		 *
		 * @param s The string to print to the profile.
		 */
		void profile(const string & s);
	
		/**
		 * @brief Print to the profile if there is one and end line.
		 *
		 * @param v The double value to print.
		 */
		void profileln(double v);
	
		/**
		 * @brief Print to the profile if there is one and end line.
		 *
		 * @param s The string to print to the profile.
		 */
		void profileln(const string & s);
	
		/**
		 * @brief Print parameters and corresponding function evaluation to profiler.
		 *
		 * @param params The parameters to print.
		 * @param value  The function evaluation.
		 */
		void printPoint(const ParameterList & params, double value);
		
		/**
		 * @brief Give a message to print to the message handler.
		 *
		 * @param message The message to print.
		 */
		void printMessage(const string & message);

		/** @} */
	
	public:
	
		static string CONSTRAINTS_AUTO;
		static string CONSTRAINTS_IGNORE;
		static string CONSTRAINTS_KEEP;
};

#endif	//_ABSTRACTOPTIMIZER_H_

