//
// File: SimpleMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: ue Nov 16 17:51 2004
//

/*
Copyright or � or Copr. CNRS, (November 19, 2004)

Julien.Dutheil@univ-montp2.fr

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

/******************************************************************************/

#include "SimpleMultiDimensions.h"

/******************************************************************************/

SimpleMultiDimensions::SimpleMultiDimensions(Function * function):
	AbstractOptimizer(function)
{
	_defaultStopCondition = new FunctionStopCondition(this);
	_stopCondition = _defaultStopCondition;
	_nbParams = 0;
	_optimizer = new BrentOneDimension(function);
}

/******************************************************************************/

SimpleMultiDimensions::~SimpleMultiDimensions()
{
	//for(unsigned int i = 0; i < _nbParams; i++) {
	//	delete _optimizers[i];
	//}
	delete _optimizer;
	delete _defaultStopCondition;
}

/******************************************************************************/

void SimpleMultiDimensions::setFunction(Function * function)
{
	AbstractOptimizer::setFunction(function);
	_optimizer -> setFunction(function);
	_stopCondition -> init();
}

/******************************************************************************/

void SimpleMultiDimensions::init(const ParameterList & params) throw (Exception)
{
	// Some cleaning first.
	// This is useful only if the optimizer have been initialized once before this time.
	//for(unsigned int i = 0; i < _nbParams; i++) {
	//	delete _optimizers[i];
	//}
	

	_parameters = params;

	_nbParams = params.size();
	if(_nbParams == 0) return;

	// Initialize optimizers:
	unsigned int nbEvalMax = _nbEvalMax / _nbParams;
	//_optimizers.resize(_nbParams);
	//for(unsigned int i = 0; i < _nbParams; i++) {
	//	//_optimizers[i] = new BrentOneDimension(_function);
	//	_optimizers[i] = new GoldenSectionSearch(_function);
	//	_optimizers[i] -> setMaximumNumberOfEvaluations(nbEvalMax);
	//	_optimizers[i] -> setProfiler(_profiler);
	//	_optimizers[i] -> setMessageHandler(_messageHandler);
	//	//_optimizers[i] -> setTolerance(_stopCondition -> getTolerance());
	//	_optimizers[i] -> getStopCondition() -> setTolerance(getStopCondition() -> getTolerance());
	//	_optimizers[i] -> setConstraintPolicy(_constraintPolicy);
	//	_optimizers[i] -> setInitialInterval(0.,1.);
	//	profile(_parameters[i] -> getName() + "\t"); 
	//}
	_optimizer -> setMaximumNumberOfEvaluations(nbEvalMax);
	_optimizer -> setProfiler(_profiler);
	_optimizer -> setMessageHandler(_messageHandler);
	_optimizer -> getStopCondition() -> setTolerance(getStopCondition() -> getTolerance());
	_optimizer -> setConstraintPolicy(_constraintPolicy);
	_optimizer -> setInitialInterval(0.,1.);
	
	for(unsigned int i = 0; i < _nbParams; i++) {
		profile(_parameters[i] -> getName() + "\t"); 
	}
	profileln("Function");

	printPoint(_parameters, _function -> f(_parameters));
	// Initialize the StopCondition:
	_stopCondition -> init();
}

/******************************************************************************/

double SimpleMultiDimensions::step() throw (Exception)
{
	double f = _function -> getValue();
	for(unsigned int i = 0; i < _nbParams; i++) {
		if(_verbose > 0) {
			cout << _parameters[i] -> getName() << ":";
			cout.flush();
		}
		// Re-init optimizer according to new values:
		double v = _parameters[i] -> getValue();
		//_optimizers[i] -> setInitialInterval(v - 0.01, v + 0.01);
		//_optimizers[i] -> init(_parameters.subList(i));
		_optimizer -> setVerbose( _verbose > 0 ? _verbose - 1 : 0);
		_optimizer -> setInitialInterval(v - 0.01, v + 0.01);
		_optimizer -> init(_parameters.subList(i));

		// Optimize through this dimension:
		//_optimizers[i] -> optimize();
		f = _optimizer -> optimize();
		// Update parameters with the new value:
		//_parameters.setParametersValues(_optimizers[i] -> getParameters());
		//_nbEval += _optimizers[i] -> getNumberOfEvaluations(); 
		_parameters.setParametersValues(_optimizer -> getParameters());
		_nbEval += _optimizer -> getNumberOfEvaluations(); 
	}
	_tolIsReached = _nbParams <= 1 || _stopCondition -> isToleranceReached();
	return f;
}

/******************************************************************************/

double SimpleMultiDimensions::optimize() throw (Exception)
{
	_tolIsReached = false;
	_nbEval = 0;
	for (_nbEval = 0; _nbEval < _nbEvalMax && ! _tolIsReached; _nbEval++) {
		step();
	}
	return _function -> getValue();
}

/******************************************************************************/

