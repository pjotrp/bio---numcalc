//
// File: Optimizer.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov  4 16:01:27 2003
//

#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "Functions.h"
#include "ParameterList.h"
#include "OptimizationStopCondition.h"

// From the STL:
#include <iostream>

using namespace std;

/**
 * @brief This is the base interface for all optimization methods.
 * 
 * An optimizer will deal with object that implement the Optimizable interface.
 */
class Optimizer
{
	public:
		Optimizer() {}
		virtual ~Optimizer() {}
	
	public:
		
		/**
		 * @brief Set the intial values of the parameters.
		 *
		 * @param params The initial values of parameters.
		 */
		virtual void init(const ParameterList & params) throw (Exception) = 0;

		/**
		 * @brief Perform a step optimization.
		 *
		 * @return the value of the function after this step.
		 */
		virtual double step() throw (Exception) = 0;

		/**
		 * @brief Get the current parameters.
		 *
		 * @return The current parameters.
		 */
		virtual ParameterList getParameters() const = 0;
		
		/**
		 * @brief Get the corresponding function evaluation.
		 *
		 * @return The value of the function at the point specified by _parameters.
		 */
		virtual double getFunctionValue() const = 0;
		
		/**
		 * @brief Perform as many optimization steps untill the stop condition is met.
		 *
		 * @return The value of the function after optimization is completed.
		 */
		virtual double optimize() throw (Exception) = 0;
	
		/**
		 * @brief Get the whole Optimizable object.
		 *
		 * @return A const pointer toward the object to be optimized.
		 */
		virtual const Function * getFunction() const = 0;

		/**
		 * @brief Set the message handler for this optimizer.
		 *
		 * The message handler keeps all messages that the optimizer may send.
		 * The default handler is set to standard output, but you can pass any
		 * ostream object (cerr, ofstream, etc.).
		 *
		 * @param mh The message handler to use.
		 */
		virtual void setMessageHandler(ostream * mh) = 0;
		
		/**
		 * @brief Set the profiler for this optimizer.
		 *
		 * The profiler keeps all the intermediate values taken by the parameters.
		 * The default profiler is set to standard output, but you can pass any
		 * ostream object (cerr, ofstream, etc.).
		 *
		 * @param profiler The profiler to use.
		 */
		virtual void setProfiler(ostream * profiler) = 0;
		
		/**
		 * @brief Get the number of function evaluation done since call of
		 * the init function.
		 *
		 * @return the number of function evaluations.
		 */
		virtual	int getNumberOfEvaluations() const = 0;
		
		/**
		 * @brief Set the stopping condition of the optimization algorithm.
		 *
		 * @param stopCondition The stoping condition to use while optimizing.
		 */
		virtual void setStopCondition(OptimizationStopCondition * stopCondition) = 0;

		/**
		 * @brief Get the stopping condition of the optimization algorithm.
		 *
		 * @return The stopping condition used while optimizing.
		 */
		virtual OptimizationStopCondition * getStopCondition() = 0;

		/**
		 * @brief Get the default stopping condition of the optimization algorithm.
		 *
		 * @return The default stopping condition used while optimizing.
		 */
		virtual OptimizationStopCondition * getDefaultStopCondition() = 0;
		
		/**
		 * @brief Tell if the tolerance level is reached.
		 *
		 * @return Whether the tolerance is reached or not.
		 * @see OptimizationStopCondition
		 */
		virtual bool isToleranceReached() const = 0;
		
		/**
		 * @brief Tell if the maximum number of function evaluations is reached.
		 *
		 * @return Whether the maximum number of function evaluations is reached or not.
		 */
		virtual bool isMaximumNumberOfEvaluationsReached() const = 0;

};

#endif	//_OPTIMIZER_H_
