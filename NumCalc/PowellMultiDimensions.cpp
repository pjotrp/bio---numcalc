//
// File: PowellMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 17 15:16:45 2003
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

#include "PowellMultiDimensions.h"

#include "NumTools.h"
#include "BrentOneDimension.h"
#include "OneDimensionOptimizationTools.h"

using namespace NumTools;

/******************************************************************************/

bool PowellMultiDimensions::PMDStopCondition::isToleranceReached() const
{
  // NRC Test for done:
  const PowellMultiDimensions * pmd = dynamic_cast<const PowellMultiDimensions *>(_optimizer);
  double fp   = pmd->_fp;
  double fret = pmd->_fret;
  return 2.0 * NumTools::abs(fp - fret) <= _tolerance * (NumTools::abs(fp) + NumTools::abs(fret));
}
    
/******************************************************************************/

PowellMultiDimensions::PowellMultiDimensions(Function * function) :
AbstractOptimizer(function), _f1dim(function)
{
  _defaultStopCondition = new PMDStopCondition(this);
  _stopCondition = _defaultStopCondition->clone();
}

/******************************************************************************/

void PowellMultiDimensions::init(const ParameterList & params) throw (Exception)
{
  AbstractOptimizer::init(params);
  
  // Build the initial matrix:
  unsigned int n = params.size();
  _xi.resize(n);
  for(unsigned int i = 0; i < n; i++)
  {
    // Copy the parameter list:
    _xi[i] = Vdouble(n);
    for(unsigned int j = 0; j < n; j++)
    {
      // Set the directions to unit vectors:
      _xi[i][j] = (j == i) ? 1 : 0;
    }
  }
  
  // Starting point:
  _fret = _function->f(_parameters);
  _pt   = _parameters;

  for (unsigned int j = 0; j < n; j++)
  {
    profile(_parameters[j]->getName() + "\t"); 
  }
  profileln("Function\tTime");
  printPoint(_parameters, _fret);
}
  
/******************************************************************************/
  
double PowellMultiDimensions::step() throw (Exception)
{
  if(!_isInitialized) throw Exception("PowellMultiDimensions::step. Optimizer not initialized: call the 'init' method first!");
  if(_verbose > 0) { cout << "*" << endl; }
  
  unsigned int n = _parameters.size();
  _fp = _fret;
  unsigned int ibig = 0;
  double del = 0.0; // Will be the biggest function decrease
  Vdouble xit(n);
  
  // In each iteration, loop over all directions in the set.
  double fptt;
  for(unsigned int i = 0; i < n; i++)
  {
    // Copy the direction:
    for(unsigned int j = 0; j < n; j++)
    {
      xit[j] = _xi[j][i];
      //xit[j] = _xi[i][j];
    }
    fptt = _fret;
    _nbEval += OneDimensionOptimizationTools::lineMinimization(_f1dim, _parameters, xit, _stopCondition->getTolerance(), _profiler, _messageHandler, _verbose);
    //_fret = _f1dim.getValue();
    _fret = _function->f(_parameters);
    printPoint(_parameters, _fret);
    if (fptt - _fret > del)
    {
      del = fptt - _fret;
      ibig = i;
    }
  }

  // Construct the extrapolated point and the average direction moved.
  // Save the old starting point.
  ParameterList ptt = _parameters;
  for(unsigned int j = 0; j < n; j++)
  {
    ptt[j]->setValue(2.0 * _parameters[j]->getValue() - _pt[j]->getValue());
    xit[j] = _parameters[j]->getValue() - _pt[j]->getValue();
    _pt[j]->setValue(_parameters[j]->getValue());
  }
  // Function value at extrapolated point.
  fptt = _function->f(ptt);
  if (fptt < _fp)
  {
    double t = 2.0 * (_fp - 2.0 * _fret + fptt) * sqr(_fp - _fret - del) - del * sqr(_fp - fptt);
    if (t < 0.0)
    {
      //cout << endl << "New direction: drection " << ibig << " removed." << endl;
      // Move to the minimum of the new direction, and save the new direction.
      _nbEval += OneDimensionOptimizationTools::lineMinimization(_f1dim, _parameters, xit, _stopCondition->getTolerance(), _profiler, _messageHandler, _verbose);
      for(unsigned int j = 0; j < n; j++)
      {
        _xi[j][ibig]  = _xi[j][n - 1];
        _xi[j][n - 1] = xit[j];
      }
    }
  }

  // We check for tolerance only once all directions have been looped over:
  _tolIsReached = _stopCondition->isToleranceReached();

  _stopCondition->init();
  return _fret;
}

/******************************************************************************/

double PowellMultiDimensions::optimize() throw (Exception)
{
  _tolIsReached = false;
  for (_nbEval = 0; _nbEval < _nbEvalMax && ! _tolIsReached; _nbEval++)
  {
    step();
  }
  // Apply best parameter:
  return _function->f(_parameters);
}

/******************************************************************************/

double PowellMultiDimensions::getFunctionValue() const throw (NullPointerException)
{
    if(_function == NULL) throw NullPointerException("PowellMultiDimensions::getFunctionValue. No function associated to this optimizer.");
    return _fret;
}

/******************************************************************************/

