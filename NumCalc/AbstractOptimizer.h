//
// File: AbstractOptimizer.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Dec 22 12:18:09 2003
//

#ifndef _ABSTRACTOPTIMIZER_H_
#define _ABSTRACTOPTIMIZER_H_

#include "Optimizer.h"

/**
 * @brief This class offers a basic implementation of the Optimizer interface.
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
		 * <ul>
		 * <li>CONSTRAINTS_KEEP: keep the constraint associated to the parameters (default).</li>
		 * <li>CONSTRAINTS_IGNORE: remove all constraints.</li>
		 * <li>CONSTRAINTS_AUTO: use AutoParameters to deal with constraints.
		 * </ul>
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


	public:
		AbstractOptimizer(Function * function);
		virtual ~AbstractOptimizer();
	
	public:
		
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */
		void init(const ParameterList & params) throw (Exception);
		ParameterList getParameters() const;
		const Function * getFunction() const;
		void setMessageHandler(ostream * mh);
		void setProfiler(ostream * profiler);
		int getNumberOfEvaluations() const;
		void setStopCondition(OptimizationStopCondition * stopCondition);
		OptimizationStopCondition * getStopCondition();
		OptimizationStopCondition * getDefaultStopCondition();
		bool isToleranceReached() const;
		bool isMaximumNumberOfEvaluationsReached() const;
		/** @} */
	
	protected: // Some util:
		
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
	
	public:

		/**
		 * @brief Set the constraint policy for this optimizer.
		 *
		 * @param constraintPolicy The constraint policy.
		 */
		void setConstraintPolicy(const string & constraintPolicy);

		/**
		 * @brief Get the constraint policy for this optimizer.
		 *
		 * @return The constraint policy.
		 */
		string getConstraintPolicy() const;
		
		/**
		 * @brief Set the maximum number of function evaluation to perform during optimization.
		 *
		 * @param max The maximum number of evaluations to perform.
		 */
		void setMaximumNumberOfEvaluations(int max);
	
	public:
	
		static string CONSTRAINTS_AUTO;
		static string CONSTRAINTS_IGNORE;
		static string CONSTRAINTS_KEEP;
};

#endif	//_ABSTRACTOPTIMIZER_H_
