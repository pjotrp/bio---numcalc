//
// File: AbstractOptimizer.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Dec 22 12:18:09 2003
//

#include "AbstractOptimizer.h"
#include "AutoParameter.h"

// From the STL:
#include <iomanip>

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

AbstractOptimizer::AbstractOptimizer(Function * function): _function(function) 
{
	// Initialization with defaults:
	_messageHandler   = & cout;
	_profiler         = & cout;
	_constraintPolicy = CONSTRAINTS_KEEP;	
	_nbEvalMax = 1000000;
}

/******************************************************************************/

AbstractOptimizer::~AbstractOptimizer() {}
	
/******************************************************************************/
	
void AbstractOptimizer::init(const ParameterList & params) throw (Exception) {
	_parameters = params;
	     if(_constraintPolicy == CONSTRAINTS_AUTO)   autoParameter();
	else if(_constraintPolicy == CONSTRAINTS_IGNORE) ignoreConstraints();
	_tolIsReached = false;
}

/******************************************************************************/

ParameterList AbstractOptimizer::getParameters() const { return _parameters; }

/******************************************************************************/

const Function * AbstractOptimizer::getFunction() const { return _function; }

/******************************************************************************/

void AbstractOptimizer::setMessageHandler(ostream * mh) { _messageHandler = mh; }

/******************************************************************************/

void AbstractOptimizer::setProfiler(ostream * profiler) { _profiler = profiler; }

/******************************************************************************/

void AbstractOptimizer::setStopCondition(
	OptimizationStopCondition * stopCondition)
{
	_stopCondition = stopCondition;
}

/******************************************************************************/

OptimizationStopCondition * AbstractOptimizer::getStopCondition() {
	return _stopCondition;
}

/******************************************************************************/

OptimizationStopCondition * AbstractOptimizer::getDefaultStopCondition() {
	return _defaultStopCondition;
}

/******************************************************************************/

bool AbstractOptimizer::isToleranceReached() const {
	return _tolIsReached;
}

/******************************************************************************/

bool AbstractOptimizer::isMaximumNumberOfEvaluationsReached() const {
	return _nbEval >= _nbEvalMax;
}

/******************************************************************************/

void AbstractOptimizer::profile(double v) {
	if(_profiler != NULL) (* _profiler) << v;
}	

/******************************************************************************/

void AbstractOptimizer::profileln(double v) {
	if(_profiler != NULL) (* _profiler) << v << endl;
}
	
/******************************************************************************/

void AbstractOptimizer::profile(const string & s) {
	if(_profiler != NULL) (* _profiler) << s;
}	

/******************************************************************************/

void AbstractOptimizer::profileln(const string & s) {
	if(_profiler != NULL) (* _profiler) << s << endl;
}
	
/******************************************************************************/

void AbstractOptimizer::printPoint(const ParameterList & params, double value) {
	int ndim = params.size();
	for (int j = 0; j < ndim; j++) {
		profile(TextTools::toString(params[j] -> getValue()));
		profile("\t"); 
	}
	profileln(value);
}

/******************************************************************************/

void AbstractOptimizer::printMessage(const string & message) {
	if(_messageHandler != NULL) (* _messageHandler) << message << endl;
}

/******************************************************************************/

void AbstractOptimizer::autoParameter() {
	for(unsigned int i = 0; i < _parameters.size(); i++) {
		Parameter * p = _parameters[i];
		AutoParameter * ap = new AutoParameter(* p);
		ap -> setMessageHandler(_messageHandler);
		_parameters[i] = ap;
		delete p;
	}
}

/******************************************************************************/

void AbstractOptimizer::ignoreConstraints() {
	for(unsigned int i = 0; i < _parameters.size(); i++) {
		_parameters[i] -> removeConstraint();
	}
}

/** Constraint policy: ********************************************************/

void AbstractOptimizer::setConstraintPolicy(const string & constraintPolicy)
{ _constraintPolicy = constraintPolicy; }

string AbstractOptimizer::getConstraintPolicy() const { return _constraintPolicy; }

string AbstractOptimizer::CONSTRAINTS_AUTO   = "auto";
string AbstractOptimizer::CONSTRAINTS_IGNORE = "ignore";
string AbstractOptimizer::CONSTRAINTS_KEEP   = "keep";

/******************************************************************************/

void AbstractOptimizer::setMaximumNumberOfEvaluations(int max) { _nbEvalMax = max; }

/******************************************************************************/

int AbstractOptimizer::getNumberOfEvaluations() const { return _nbEval; }

/******************************************************************************/
