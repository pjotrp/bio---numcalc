//
// File: GoldenSectionSearch.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 10 10:42:17 2003
//

#include "GoldenSectionSearch.h"
#include "NumTools.h"
#include "OneDimensionOptimizationTools.h"

using namespace NumTools;

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

GoldenSectionSearch::GSSStopCondition::GSSStopCondition(GoldenSectionSearch * gss):
	AbstractOptimizationStopCondition(gss) {}

GoldenSectionSearch::GSSStopCondition::~GSSStopCondition() {}
		
bool GoldenSectionSearch::GSSStopCondition::isToleranceReached() const {
	// NRC Test for done:
	const GoldenSectionSearch * gss = dynamic_cast<const GoldenSectionSearch *>(_optimizer);
	double x0 = gss -> x0;
	double x1 = gss -> x1;
	double x2 = gss -> x2;
	double x3 = gss -> x3;
	return NumTools::abs(x3 - x0) > _tolerance * (NumTools::abs(x1) + NumTools::abs(x2));
}
		
/******************************************************************************/

double GoldenSectionSearch::R = 0.61803399;
double GoldenSectionSearch::C = 1. - R;

/******************************************************************************/

GoldenSectionSearch::GoldenSectionSearch(const Function * function) :
AbstractOptimizer(function)
{
	_nbEvalMax = 10000;
	_defaultStopCondition = new GSSStopCondition(this);
	_stopCondition = _defaultStopCondition;
}

/******************************************************************************/

GoldenSectionSearch::~GoldenSectionSearch() {
	delete _defaultStopCondition;
}

/******************************************************************************/

void GoldenSectionSearch::init(const ParameterList & params) throw (Exception) {
	// Set the initial value (no use here! Use setInitialValues() instead).
	if(params.size() != 1) throw Exception("GoldenSectionSearch::init(). This optimizer only deals with one parameter.");
	AbstractOptimizer::init(params);

	// Bracket the minimum.
	Bracket bracket = OneDimensionOptimizationTools::bracketMinimum(_xinf, _xsup, _function, _parameters);
	//printMessage("Initial bracketing:");
	//printMessage("A: x = " + TextTools::toString(bracket.a.x) + ", f = " + TextTools::toString(bracket.a.f));
	//printMessage("B: x = " + TextTools::toString(bracket.b.x) + ", f = " + TextTools::toString(bracket.b.f));
	//printMessage("C: x = " + TextTools::toString(bracket.c.x) + ", f = " + TextTools::toString(bracket.c.f));

	// At any given time we will keep track of four points, x0, x1, x2 and x3.
	x0 = bracket.a.x;
	x3 = bracket.c.x;
	if (NumTools::abs(bracket.c.x - bracket.b.x)
      > NumTools::abs(bracket.b.x - bracket.a.x)) {
		// Make x0 to x1 the smaller segment,
		x1 = bracket.b.x;
		// and fill in the new point to be tried.
		x2 = bracket.b.x + C * (bracket.c.x - bracket.b.x);
	} else {
		x2 = bracket.b.x;
		x1 = bracket.b.x - C * (bracket.b.x - bracket.a.x);
	}
	// The initial function evaluations.
	// Note that we never need to evaluate the function at the original endpoints.
	_parameters[0] -> setValue(x1); f1 = _function -> f(_parameters);
	_parameters[0] -> setValue(x2); f2 = _function -> f(_parameters);
	_nbEval = 0;	

	profileln(_parameters[0] -> getName() + "\tFunction");
}

/******************************************************************************/

void GoldenSectionSearch::setInitialInterval(double inf, double sup) {
	_xinf = inf; _xsup = sup;	
}

/******************************************************************************/

double GoldenSectionSearch::step() throw (Exception)
{
	cout << "."; cout.flush();
	
	_nbEval++;

	if (f2 < f1) {
		// One possible outcome, its housekeeping,
		NumTools::shift<double>(x0, x1, x2);
		x2 = R * x1 + C * x3;
		// and a new function evaluation.
		_parameters[0] -> setValue(x2);
		_tolIsReached = _nbEval > 2 && _stopCondition -> isToleranceReached();
		NumTools::shift<double>(f1, f2, _function -> f(_parameters));
		printPoint(_parameters, f2);
		return f2;
	} else {
		// The other outcome,
		NumTools::shift<double>(x3, x2, x1);
		x1 = R * x2 + C * x0;
		// and its new function evaluation.
		_parameters[0] -> setValue(x1);
		_tolIsReached = _stopCondition -> isToleranceReached();
		NumTools::shift<double>(f2, f1, _function -> f(_parameters));
		printPoint(_parameters, f1);
		return f1;
	}
}

/******************************************************************************/

double GoldenSectionSearch::optimize() throw (Exception) {
	_tolIsReached = false;
	while (_nbEval < _nbEvalMax && !_tolIsReached) {
		step();
	} // Back to see if we are done.

	// We are done. Output the best of the two current values.
	// Apply the min value to the function:
	if (f1 < f2) {
		_parameters[0] -> setValue(x1); f1 = _function -> f(_parameters);
		return f1;
	} else {
		_parameters[0] -> setValue(x2); f2 = _function -> f(_parameters);
		return f2;
	}
}

/******************************************************************************/

double GoldenSectionSearch::getFunctionValue() const {
	return NumTools::min(f1, f2); 
}

/******************************************************************************/
