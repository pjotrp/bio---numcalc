//
// File: DownhillSimplexMethod.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov  4 17:10:05 2003
//

#include "DownhillSimplexMethod.h"
#include "NumTools.h"
using namespace NumTools;

/******************************************************************************/

DownhillSimplexMethod::DSMStopCondition::DSMStopCondition(DownhillSimplexMethod * dsm):
	AbstractOptimizationStopCondition(dsm) {}
		
DownhillSimplexMethod::DSMStopCondition::~DSMStopCondition() {}

bool DownhillSimplexMethod::DSMStopCondition::isToleranceReached() const
{
	// NRC stop condition, replaced by a general stop condition on parmeter estimates.
	const DownhillSimplexMethod * dsm = dynamic_cast<const DownhillSimplexMethod *>(_optimizer);
	Vdouble y    = dsm -> _y;
	int iLowest  = dsm -> _iLowest;
	int iHighest = dsm -> _iHighest;
	// Compute the fractional range from highest to lowest and return if satisfactory.
	double rTol = 2.0 * NumTools::abs(y[iHighest] - y[iLowest]) /
		(NumTools::abs(y[iHighest]) + NumTools::abs(y[iLowest]));
	return rTol < _tolerance;
}
	
/******************************************************************************/

DownhillSimplexMethod::Simplex::Simplex(): vector<ParameterList>() {}
DownhillSimplexMethod::Simplex::~Simplex() {}
			
int DownhillSimplexMethod::Simplex::getDimension() const { return operator[](0).size(); }

/******************************************************************************/
			
DownhillSimplexMethod::DownhillSimplexMethod(const Function * function): AbstractOptimizer(function) {
	// Default values:
	_nbEvalMax = 5000;
	_defaultStopCondition = new DSMStopCondition(this);
	_stopCondition = _defaultStopCondition;
}

/******************************************************************************/
	
DownhillSimplexMethod::~DownhillSimplexMethod() {
		delete _defaultStopCondition;
}

/******************************************************************************/

void DownhillSimplexMethod::init(const ParameterList & params) throw (Exception)
{
	AbstractOptimizer::init(params);

	int nDim = _parameters.size();

	// Initialize the simplex:
	_simplex.resize(nDim + 1);
	_y.resize(nDim + 1);
	_simplex[0] = _parameters;
	_y[0] = _function -> f(_simplex[0]);
	double lambda = 1.;
	for(int i = 1; i < nDim + 1; i++) {
		// Copy the vector...
		_simplex[i] = _parameters;
		// ... and set the initial values.
		for(int j = 0; j < nDim; j++) {
			_simplex[i][j] -> setValue(_parameters[j] -> getValue() + (j == i - 1 ? lambda : 0.));
		}
		//Compute the corresponding f value:
		_y[i] = _function -> f(_simplex[i]);
	}
	
	_nbEval = 0;
	_pSum = getPSum();

	for (int j = 0; j < nDim; j++) {
		profile(_parameters[j] -> getName() + "\t"); 
	}
	profileln("Function");
}
	
/******************************************************************************/

double DownhillSimplexMethod::step() throw (Exception) {
	cout << "."; cout.flush();
	// The number of dimensions of the parameter space:
	int nDim = _simplex.getDimension();
	int mpts = nDim + 1;

	_iLowest = 0;
	// First we must determine which point is the highest (worst),
	// next-highest, and lowest (best), by looping over the points
	// in the simplex.
	if (_y[0] > _y[1]) {
		_iHighest = 0;
		_iNextHighest = 1;
	} else {
		_iHighest = 1;
		_iNextHighest = 0;
	}
	
	for (int i = 0; i < mpts; i++) {
		if (_y[i] <= _y[_iLowest]) _iLowest = i;
		if (_y[i] > _y[_iHighest]) {
			_iNextHighest = _iHighest;
			_iHighest = i;
		} else if (_y[i] > _y[_iNextHighest] && i != _iHighest) _iNextHighest = i;
	}
		
	_nbEval += 2;

	// Begin a new iteration.
	// First extrapolate by a factor -1 through the face of the simplex
	// across from high point, i.e., reflect the simplex from the high point.</p>

	double yTry = amotry(-1.0);
	if (yTry <= _y[_iLowest]) {
		// Gives a result better then the best point,
		// so try an additional extrapolation by a factor 2.
		yTry = amotry(2.0);

		// Set current best point:
		_parameters = _simplex[_iLowest];

		// Test for stop:
		_tolIsReached = _nbEval > 2 && _stopCondition -> isToleranceReached();

		// Print parameters to profile:
		printPoint(_simplex[_iLowest], _y[_iLowest]);
	} else if (yTry >= _y[_iNextHighest]) {
		// The reflect point is worse than the second-highest,
		// so look for an intermediate lower point, i.e., do a one-dimensional
		// contraction.
		double ySave = _y[_iHighest];
		yTry = amotry(0.5);
		if (yTry >= ySave) {
			// Can't seem to get rid of that high point.
			// Better contract around the lowest (best) point.
			for (int i = 0; i < mpts; i++) {
				if (i != _iLowest) {
					for (int j = 0; j < nDim; j++) {
						_pSum[j] -> setValue(0.5 * (_simplex[i][j] -> getValue() + _simplex[_iLowest][j] -> getValue()));
						_simplex[i][j] -> setValue(_pSum[j] -> getValue());
					}
					_y[i] = _function -> f(_pSum);
					//printPoint(_pSum, _y[i]);
				}
			}
			_nbEval += nDim;
			_pSum = getPSum();
		}
	} else --(_nbEval); // Correct the evaluation count.

	// Send current value of lower point:
	return _y[_iLowest];
}

/******************************************************************************/

double DownhillSimplexMethod::optimize() throw (Exception) {

	_nbEval = 0;
	_tolIsReached = false;
	while (_nbEval < _nbEvalMax && !_tolIsReached) {
		step();
	}

	// set best shot:
	return _function -> f(_simplex[_iLowest]);
}

/******************************************************************************/

double DownhillSimplexMethod::getFunctionValue() const { return _y[_iLowest]; }

/******************************************************************************/

ParameterList DownhillSimplexMethod::getPSum() {
	int ndim = _simplex.getDimension();
	int mpts = ndim + 1;
	
	// Get a copy...
	ParameterList pSum = _parameters;
	// ... and initializes it.
	for (int j = 0; j < ndim; j++) {
		double sum = 0.;
		for (int i = 0; i < mpts; i++) {
			sum += _simplex[i][j] -> getValue();
		}
		pSum[j] -> setValue(sum);
	}
	return pSum;
}

/******************************************************************************/

double DownhillSimplexMethod::amotry(double fac)
{
	int ndim = _simplex.getDimension();
	double fac1, fac2, yTry;

	fac1 = (1.0 - fac) / ndim;
	fac2 = fac1 - fac;
	
	// Get a copy...
	ParameterList pTry = _parameters;
	// and initialize it:
	for (int j = 0; j < ndim; j++) {
		pTry[j] -> setValue(_pSum[j] -> getValue() * fac1 - _simplex[_iHighest][j] -> getValue() * fac2);
	}
	// Now compute the function for this new set of parameters:
	yTry = _function -> f(pTry);
	
	// Then test this new point:
	if (yTry < _y[_iHighest]) {
		_y[_iHighest] = yTry;
		for (int j = 0; j < ndim; j++) {
			_pSum[j] -> setValue(_pSum[j] -> getValue() + pTry[j] -> getValue() - _simplex[_iHighest][j] -> getValue());
			_simplex[_iHighest][j] -> setValue(pTry[j] -> getValue());
		}
	}
	return yTry;
}

/******************************************************************************/
