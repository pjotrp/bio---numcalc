//
// File: OptimizationStopCondition.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Dec 23 11:51:31 2003
//

#include "OptimizationStopCondition.h"
#include "Optimizer.h"
#include "VectorTools.h"
using namespace VectorFunctions;
#include "NumTools.h"
using namespace NumTools;

/******************************************************************************/

OptimizationStopCondition::OptimizationStopCondition() {}

OptimizationStopCondition::~OptimizationStopCondition() {}

/******************************************************************************/

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(const Optimizer * optimizer):
	_optimizer(optimizer),
	_tolerance(0.),
	_callCount(0),
	_burnin(0) {}

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(
	const Optimizer * optimizer,
	double tolerance):
	_optimizer(optimizer),
	_tolerance(tolerance),
	_callCount(0),
	_burnin(0) {}

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(
	const Optimizer * optimizer,
	int burnin):
	_optimizer(optimizer),
	_tolerance(0.),
	_callCount(0),
	_burnin(burnin) {}

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(
	const Optimizer * optimizer,
	double tolerance,
	int burnin):
	_optimizer(optimizer),
	_tolerance(tolerance),
	_callCount(0),
	_burnin(burnin) {}

AbstractOptimizationStopCondition::~AbstractOptimizationStopCondition() {}

void AbstractOptimizationStopCondition::setTolerance(double tolerance) {
	_tolerance = tolerance;
}

double AbstractOptimizationStopCondition::getTolerance() const {
	return _tolerance;
}

void AbstractOptimizationStopCondition::resetCounter() {
	_callCount = 0;
}

void AbstractOptimizationStopCondition::setBurnin(int burnin) {
	_burnin = burnin;
}

int AbstractOptimizationStopCondition::getBurnin() const {
	return _burnin;
}

/******************************************************************************/

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer):
	AbstractOptimizationStopCondition(optimizer)
{
	_newParametersEstimates = _optimizer -> getParameters();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer,
	double tolerance):
	AbstractOptimizationStopCondition(optimizer, tolerance)
{
	_newParametersEstimates = _optimizer -> getParameters();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer,
	int burnin):
	AbstractOptimizationStopCondition(optimizer, burnin)
{
	_newParametersEstimates = _optimizer -> getParameters();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer,
	double tolerance,
	int burnin):
	AbstractOptimizationStopCondition(optimizer, tolerance, burnin)
{
	_newParametersEstimates = _optimizer -> getParameters();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::~ParametersStopCondition() {}

/******************************************************************************/

bool ParametersStopCondition::isToleranceReached() const
{
	_callCount++;
	_lastParametersEstimates = _newParametersEstimates;
	_newParametersEstimates   = _optimizer -> getParameters();
	if(_callCount <= _burnin) return false;
	for(unsigned int i = 0; i < _newParametersEstimates.size(); i++) {
		Parameter * p = _newParametersEstimates[i];
		if(p == NULL) throw ParameterNotFoundException("ParameterStopCondition::isToleranceReached.", p -> getName());
		double lastEstimate = _lastParametersEstimates.getParameter(p -> getName()) -> getValue();
		double newEstimate = p -> getValue();
		double tol = NumTools::abs<double>(newEstimate - lastEstimate);
		if(tol > _tolerance) {
			return false;
		}
	}
	return true;
}

/******************************************************************************/
