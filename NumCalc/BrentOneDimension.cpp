//
// File: BrentOneDimension.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 17 11:45:58 2003
//

#include "BrentOneDimension.h"
#include "NumTools.h"
#include "OneDimensionOptimizationTools.h"

using namespace NumTools;

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

BrentOneDimension::BODStopCondition::BODStopCondition(BrentOneDimension * bod):
	AbstractOptimizationStopCondition(bod) {}

BrentOneDimension::BODStopCondition::~BODStopCondition() {}
		
bool BrentOneDimension::BODStopCondition::isToleranceReached() const {
	// NRC Test for done:
	const BrentOneDimension * bod = dynamic_cast<const BrentOneDimension *>(_optimizer);
	double x    = bod -> x;
	double xm   = bod -> xm;
	double tol2 = bod -> tol2;
	double b    = bod -> b;
	double a    = bod -> a;
	return NumTools::abs(x - xm) <= (tol2 - 0.5 * (b - a));
}
		
/******************************************************************************/

BrentOneDimension::BrentOneDimension(const Function * function) :
AbstractOptimizer(function)
{
	_defaultStopCondition = new BODStopCondition(this);
	_stopCondition = _defaultStopCondition;
	_nbEvalMax = 10000;
}


BrentOneDimension::~BrentOneDimension()
{
	delete _defaultStopCondition;
}
	
/******************************************************************************/
	
double BrentOneDimension::CGOLD = 0.3819660;
double BrentOneDimension::ZEPS  = 1.0e-10;
	
/******************************************************************************/
	
void BrentOneDimension::init(const ParameterList & params) throw (Exception) {
	// Set the initial value (no use here! Use setInitialValues() instead).
	if(params.size() != 1) throw Exception("BrentOneDimension::init(). This optimizer only deals with one parameter.");
	AbstractOptimizer::init(params);

	// Bracket the minimum.
	Bracket bracket = OneDimensionOptimizationTools::bracketMinimum(_xinf, _xsup, _function, _parameters);
	//printMessage("Initial bracketing:");
	//printMessage("A: x = " + TextTools::toString(bracket.a.x) + ", f = " + TextTools::toString(bracket.a.f));
	//printMessage("B: x = " + TextTools::toString(bracket.b.x) + ", f = " + TextTools::toString(bracket.b.f));
	//printMessage("C: x = " + TextTools::toString(bracket.c.x) + ", f = " + TextTools::toString(bracket.c.f));

	// This will be the distance moved on the step before last.
	e = 0.0;

	// a and b must be in ascending order, but input abscissa need not be.
	a = (bracket.a.x < bracket.c.x ? bracket.a.x : bracket.c.x);
	b = (bracket.a.x > bracket.c.x ? bracket.a.x : bracket.c.x);
	
	// Initializations...
	x = w = v = bracket.b.x;
	_parameters[0] -> setValue(x);
	fw = fv = fx = _function -> f(_parameters);	

	profileln(_parameters[0] -> getName() + "\tFunction");
}

/******************************************************************************/

void BrentOneDimension::setInitialInterval(double inf, double sup) {
	_xinf = inf; _xsup = sup;	
}

/******************************************************************************/

void BrentOneDimension::setTolerance(double tolerance) {
	_tolerance = tolerance;
}

/******************************************************************************/

double BrentOneDimension::step() throw (Exception) {
	cout << "."; cout.flush();
	xm   = 0.5 * (a + b);
	tol2 = 2.0 * (tol1 = _tolerance * NumTools::abs(x) + ZEPS);
	
	// Construct a trial parabolic fit:
	if (NumTools::abs(e) > tol1) {
		r = (x - w) * (fx - fv);
		q = (x - v) * (fx - fw);
		p = (x - v) * q - (x - w) * r;
		q = 2.0 * (q - r);
		if (q > 0.0) p = -p;
		q = NumTools::abs(q);
		etemp = e;
		e = d;
		if (NumTools::abs(p) >= NumTools::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		// The above conditions determine the acceptability of the parabolic fit.
		// Here we take the golden section step into the larger of the two segments.
		else {
			d = p / q; // Take the parabolic step.
			u = x + d;
			if (u - a < tol2 || b - u < tol2)
				d = sign(tol1, xm - x);
		}
	} else {
		d = CGOLD * (e = (x >= xm ? a - x : b - x));
	}
	u = (NumTools::abs(d) >= tol1 ? x + d : x + sign(tol1, d));

	// This is the one function evaluate per iteration.
	ParameterList pl = _parameters;
	pl[0] -> setValue(u);
	fu = _function -> f(pl);
	// Now decide what to do with our function evaluation.
	if (fu <= fx) {
		if (u >= x) a = x; else b = x;
		// Here is the house keeping:
		shift(v, w, x, u);
		shift(fv, fw, fx, fu);
	} else {
		if (u < x) a = u; else b = u;
		if (fu <= fw || w == x) {
			v  = w;
			w  = u;
			fv = fw;
			fw = fu;
		} else if (fu <= fv || v == x || v == w) {
			v  = u;
			fv = fu;
		}
	} // Done with housekeeping.

	// Store results for this step:
	_parameters[0] -> setValue(x);
	_tolIsReached = _nbEval > 2 && _stopCondition -> isToleranceReached();
	printPoint(_parameters, fx);
	return fx;
}

/******************************************************************************/
	
double BrentOneDimension::optimize() throw (Exception)
{
	// Main program loop.
	_tolIsReached = false;
	for (_nbEval = 0; _nbEval < _nbEvalMax && !_tolIsReached; _nbEval++) {
		step();
	}

	// Apply parameters and evaluate function at minimum point:
	_parameters[0] -> setValue(x); fx = _function -> f(_parameters);
	return fx;
}

/******************************************************************************/

double BrentOneDimension::getFunctionValue() const
{
	return fx;
}

/******************************************************************************/
