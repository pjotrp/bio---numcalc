//
// File: OptimizationStopCondition.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Dec 23 11:51:31 2003
//

#ifndef _OPTIMIZATIONSTOPCONDITION_H_
#define _OPTIMIZATIONSTOPCONDITION_H_

#include "ParameterList.h"

class Optimizer;
	
/******************************************************************************/

class OptimizationStopCondition
{
	public:
		OptimizationStopCondition();
		virtual ~OptimizationStopCondition();
	
	public:

		/**
		 * @brief Tell if the we reached the desired tolerance with a given 
		 * new set of estimates.
		 *
		 * The new parameter list is compared to the last estimates,
		 * and the lastParameterEstimates list is actulaized with the newParameters list.
		 *
		 * @return True if the tolerance level is reached.
		 */
		virtual bool isToleranceReached() const = 0;
		
		/**
		 * @brief Set the tolerance parameter.
		 *
		 * @param tol The tolerance parameter.
		 */	
		virtual void setTolerance(double tolerance) = 0;

		/**
		 * @brief Get the tolerance parameter.
		 *
		 * @return The tolerance parameter.
		 */	
		virtual double getTolerance() const = 0;
};

/******************************************************************************/

class AbstractOptimizationStopCondition: public OptimizationStopCondition
{
	protected:
		const Optimizer * _optimizer;
		double _tolerance;

		/**
		 * @brief Count the number of times the isToleranceReached() function
		 * has been called.
		 */
		mutable double _callCount;
	
		int _burnin;

	
	public:
		AbstractOptimizationStopCondition(const Optimizer * optimizer);
		AbstractOptimizationStopCondition(const Optimizer * optimizer, double tolerance);
		AbstractOptimizationStopCondition(const Optimizer * optimizer, int burnin);
		AbstractOptimizationStopCondition(const Optimizer * optimizer, double tolerance, int burnin);
	
		virtual ~AbstractOptimizationStopCondition();

	public:
		void setTolerance(double tolerance);
		double getTolerance() const;
	
		virtual void resetCounter();
		virtual void setBurnin(int burnin);
		virtual int getBurnin() const;

};
	
/******************************************************************************/

class ParametersStopCondition: public AbstractOptimizationStopCondition
{
	protected:
		/**
		 * @brief The last estimates of the parameters.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable ParameterList _lastParametersEstimates;
		
		/**
		 * @brief The new estimates of the parameters.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable ParameterList _newParametersEstimates;
	
	public:
		ParametersStopCondition(const Optimizer * optimizer);
		ParametersStopCondition(const Optimizer * optimizer, double tolerance);
		ParametersStopCondition(const Optimizer * optimizer, int burnin);
		ParametersStopCondition(const Optimizer * optimizer, double tolerance, int burnin);
		
		~ParametersStopCondition();
	
	public:
		bool isToleranceReached() const;
};

/******************************************************************************/

class FunctionStopCondition: public AbstractOptimizationStopCondition
{
	protected:
		/**
		 * @brief The last value of the function.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable double _lastFunctionValue;
		
		/**
		 * @brief The new value of the function.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable double _newFunctionValue;
	
	public:
		FunctionStopCondition(const Optimizer * optimizer);
		FunctionStopCondition(const Optimizer * optimizer, double tolerance);
		FunctionStopCondition(const Optimizer * optimizer, int burnin);
		FunctionStopCondition(const Optimizer * optimizer, double tolerance, int burnin);
		
		~FunctionStopCondition();
	
	public:
		bool isToleranceReached() const;

};

/******************************************************************************/

#endif	//_OPTIMIZATIONSTOPCONDITION_H_
