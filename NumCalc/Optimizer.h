//
// File: Optimizer.h
// Created by: Julien Dutheil
// Created on: Tue Nov  4 16:01:27 2003
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

#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "Functions.h"
#include "ParameterList.h"
#include "OptimizationStopCondition.h"

//From Utils:
#include <Utils/Clonable.h>

// From the STL:
#include <iostream>

using namespace std;

/**
 * @brief This is the basal interface for all optimization methods.
 * 
 * An optimizer deals with Function objects.
 */
class Optimizer: public Clonable
{
	public:
		Optimizer() {}
		virtual ~Optimizer() {}

    Optimizer * clone() const = 0;
	
	public:
		
		/**
		 * @brief Set the initial values of the parameters.
		 *
		 * @param params The initial values of parameters.
		 * @throw Exception If a problem occured during initialization.
		 */
		virtual void init(const ParameterList & params) throw (Exception) = 0;

		/**
		 * @brief Perform an optimization step.
		 *
		 * @return the value of the function after this step.
		 * @throw Exception If a problem occured during optimization.
		 */
		virtual double step() throw (Exception) = 0;

		/**
		 * @return The parameters with their current values.
		 */
		virtual ParameterList getParameters() const = 0;
		
		/**
		 * @brief Get the current function value.
		 *
		 * @return The value of the function at the point specified by _parameters.
		 * @throw NullPointerException If no function is associated with this optimizer.
		 */
		virtual double getFunctionValue() const throw (NullPointerException) = 0;
		
		/**
		 * @brief Perform as many optimization steps untill the stop condition is met.
		 *
		 * @return The value of the function after optimization is completed.
		 */
		virtual double optimize() throw (Exception) = 0;
	
		/**
		 * @brief Set the function to optimize.
		 *
		 * @param function The function to optimize.
		 */
		virtual void setFunction(Function * function) = 0;
		
		/**
		 * @brief Get the current function being optimized.
		 *
		 * @return A const pointer toward the function being optimized.
		 */
		virtual const Function * getFunction() const = 0;

		/**
		 * @brief Get the current function being optimized.
		 *
		 * @return A const pointer toward the function being optimized.
		 */
		virtual Function * getFunction() = 0;

		/**
		 * @brief Set the message handler for this optimizer.
		 *
		 * The message handler keeps all messages that the optimizer may send.
		 * The default handler is set to standard output, but you can pass any
		 * ostream object (cerr, ofstream, etc.).
		 *
		 * A NULL pointer disables message output.
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
		 * A NULL pointer disables message output.
		 * 
		 * @param profiler The profiler to use.
		 */
		virtual void setProfiler(ostream * profiler) = 0;
		
		/**
		 * @brief Get the number of function evaluations performed since the
		 * call of the init function.
		 *
		 * @return The number of function evaluations.
		 */
		virtual	int getNumberOfEvaluations() const = 0;
		
		/**
		 * @brief Set the stop condition of the optimization algorithm.
		 *
		 * @param stopCondition The stop condition to use while optimizing.
		 * @see OptimizationStopCondition.
		 */
		virtual void setStopCondition(const OptimizationStopCondition & stopCondition) = 0;

		/**
		 * @brief Get the stop condition of the optimization algorithm.
		 *
		 * @return The stop condition used while optimizing.
		 */
		virtual OptimizationStopCondition * getStopCondition() = 0;

		/**
		 * @brief Get the stop condition of the optimization algorithm.
		 *
		 * @return The stop condition used while optimizing.
		 */
		virtual const OptimizationStopCondition * getStopCondition() const = 0;

		/**
		 * @brief Get the default stop condition of the optimization algorithm.
		 *
		 * @return The default stop condition used while optimizing.
		 */
		virtual OptimizationStopCondition * getDefaultStopCondition() = 0;
		
		/**
		 * @brief Get the default stop condition of the optimization algorithm.
		 *
		 * @return The default stop condition used while optimizing.
		 */
		virtual const OptimizationStopCondition * getDefaultStopCondition() const = 0;
		
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

		/**
		 * @brief Set the maximum number of function evaluation to perform during optimization.
		 *
		 * @param max The maximum number of evaluations to perform.
		 */
		virtual void setMaximumNumberOfEvaluations(unsigned int max) = 0;

		/**
		 * @brief Set the verbose level.
		 *
		 * 0 = off
		 * 1 = on
		 * 2 = more verbose
		 * 3 = even more, etc.
		 *
		 * In most cases, only the 0 and 1 levels are implemented.
		 *
		 * @param v verbose level.
		 */
		virtual void setVerbose(unsigned int v) = 0;

		/**
		 * @brief Get the verbose level.
		 *
		 * @return verbose level.
		 */
		virtual unsigned int getVerbose() const = 0;

		/**
		 * @brief Set the constraint policy for this optimizer.
		 *
		 * @param constraintPolicy The constraint policy.
		 */
		virtual void setConstraintPolicy(const string & constraintPolicy) = 0;

		/**
		 * @brief Get the constraint policy for this optimizer.
		 *
		 * @return The constraint policy.
		 */
		virtual string getConstraintPolicy() const = 0;

};

#endif	//_OPTIMIZER_H_

