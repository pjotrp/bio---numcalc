//
// File: PowellMultiDimensions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 17 15:16:45 2003
//

#include "PowellMultiDimensions.h"

#include "NumTools.h"
#include "BrentOneDimension.h"

using namespace NumTools;

/******************************************************************************/

PowellMultiDimensions::PMDStopCondition::PMDStopCondition(PowellMultiDimensions * pmd):
	AbstractOptimizationStopCondition(pmd) {}

PowellMultiDimensions::PMDStopCondition::~PMDStopCondition() {}
		
bool PowellMultiDimensions::PMDStopCondition::isToleranceReached() const {
	// NRC Test for done:
	const PowellMultiDimensions * pmd = dynamic_cast<const PowellMultiDimensions *>(_optimizer);
	double fp   = pmd -> _fp;
	double fret = pmd -> _fret;
	return 2.0 * NumTools::abs(fp - fret) <= _tolerance * (NumTools::abs(fp) + NumTools::abs(fret));
}
		
/** Constructor: **************************************************************/

PowellMultiDimensions::PowellMultiDimensions(Function * function) :
AbstractOptimizer(function)
{
	f1dim = new PowellMultiDimensions::DirectionFunction(_function);
	_defaultStopCondition = new PMDStopCondition(this);
	_stopCondition = _defaultStopCondition;
}

/** Destructor: ***************************************************************/

PowellMultiDimensions::~PowellMultiDimensions()
{
	delete f1dim;
	delete _defaultStopCondition;
}

/******************************************************************************/

void PowellMultiDimensions::init(const ParameterList & params) throw (Exception)
{
	AbstractOptimizer::init(params);
	
	// Build the initial matrix:
	unsigned int n = params.size();
	_xi.resize(n);
	for(unsigned int i = 0; i < n; i++) {
		// Copy the parameter list:
		_xi[i] = Vdouble(n);
		for(unsigned int j = 0; j < n; j++) {
			// Set the directions to unit vectors:
			_xi[i][j] = (j == i) ? 1 : 0;
		}
	}
	
	// Starting point:
	_fret = _function -> f(_parameters);
	_pt   = _parameters;

	for (unsigned int j = 0; j < n; j++) {
		profile(_parameters[j] -> getName() + "\t"); 
	}
	profileln("Function");
	printPoint(_parameters, _fret);
}
	
/******************************************************************************/
	
double PowellMultiDimensions::step() throw (Exception)
{
	cout << "*" << endl;
	
	int n = _parameters.size();
	_fp = _fret;
	int ibig = 0;
	double del = 0.0; // Will be the biggest function decrease
	Vdouble xit(n);
	
	// In each iteration, loop over all directions in the set.
	double fptt;
	for (int i = 0; i < n; i++) {
		// Copy the direction:
		for(int j = 0; j < n; j++) {
			xit[j] = _xi[j][i];
			//xit[j] = _xi[i][j];
		}
		fptt = _fret;
		linmin(xit);
		if (fptt - _fret > del) {
			del = fptt - _fret;
			ibig = i;
		}
	}

	// Construct the extrapolated point and the average direction moved.
	// Save the old starting point.
	ParameterList ptt = _parameters;
	for (int j = 0; j < n; j++) {
		ptt[j] -> setValue(2.0 * _parameters[j] -> getValue() - _pt[j] -> getValue());
		xit[j] = _parameters[j] -> getValue() - _pt[j] -> getValue();
		_pt[j] -> setValue(_parameters[j] -> getValue());
	}
	// Function value at extrapolated point.
	fptt = _function -> f(ptt);
	if (fptt < _fp) {
		double t = 2.0 * (_fp - 2.0 * _fret + fptt) * sqr(_fp - _fret - del) - del * sqr(_fp - fptt);
		if (t < 0.0) {
			//cout << endl << "New direction: drection " << ibig << " removed." << endl;
			// Move to the minimum of the new direction, and save the new direction.
			linmin(xit);
			for (int j = 0; j < n; j++) {
				_xi[j][ibig]  = _xi[j][n - 1];
				_xi[j][n - 1] = xit[j];
			}
		}
	}

	// We check for tolerance only once all directions have been looped over:
	_tolIsReached = _stopCondition -> isToleranceReached();

	return _fret;
}

/******************************************************************************/

double PowellMultiDimensions::optimize() throw (Exception)
{
	_tolIsReached = false;
	for (_nbEval = 0; _nbEval < _nbEvalMax && ! _tolIsReached; _nbEval++) {
		step();
	}
	// Apply best parameter:
	return _function -> f(_parameters);
}

/******************************************************************************/

double PowellMultiDimensions::getFunctionValue() const { return _fret; }

/******************************************************************************/

void PowellMultiDimensions::linmin(Vdouble & xi)
{
	int n = _parameters.size();
	
	// Initial guess for brackets:
	double ax = 0.0;
	double xx = 1.0;
	
	//_parameters.printParameters(cout); cout << endl;
	f1dim -> set(_parameters, xi);
	BrentOneDimension bod(f1dim);
	bod.setMessageHandler(_messageHandler);
	bod.setProfiler(NULL);
	bod.getStopCondition() -> setTolerance(_stopCondition -> getTolerance());
	bod.setInitialInterval(ax, xx);
	bod.setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_KEEP);
	ParameterList singleParameter;
	singleParameter.addParameter(Parameter("x", 0.0));
	bod.init(singleParameter);
	_fret = bod.optimize();
	_nbEval += bod.getNumberOfEvaluations();
	
	cout << "#"; cout.flush();
	
	double xmin = f1dim -> getParameters()[0] -> getValue();
	//cout << "xmin = " << xmin << endl;
	for (int j = 0; j < n; j++) {
		//cout << "xi[" << j << "] = " << xi[j] -> getValue() * xmin << endl;
		         xi[j] *=  xmin;
		//cout << "xi[" << j << "] = " << xi[j] -> getValue() << endl;
		_parameters[j] -> setValue(_parameters[j] -> getValue() + xi[j]);
	}
	printPoint(_parameters, _fret);
}

/** DirectionFunction: ********************************************************/

PowellMultiDimensions::DirectionFunction::DirectionFunction(
	Function * function):
	function(function) {}

PowellMultiDimensions::DirectionFunction::~DirectionFunction() {}
	
void PowellMultiDimensions::DirectionFunction::setParameters(
	const ParameterList & params)
	throw (ParameterNotFoundException, ConstraintException)
{
	_params = params;
	_xt = p;
	for(unsigned int j = 0; j < p.size(); j++) {
		_xt[j] -> setValue((p[j] -> getValue()) + (_params[0] -> getValue()) * xi[j]);
	}
}

double PowellMultiDimensions::DirectionFunction::getValue() const throw (Exception)
{
	return function -> f(_xt);
}

ParameterList PowellMultiDimensions::DirectionFunction::getParameters() const throw (Exception) {
	return _params;
}

void PowellMultiDimensions::DirectionFunction::set(const ParameterList & p, const Vdouble & xi) {
	this -> p = p;
	this -> xi = xi;	
}

/******************************************************************************/
