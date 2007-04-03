//
// File: AbstractOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 22 12:18:09 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "AbstractOptimizer.h"
#include "AutoParameter.h"

// From the STL:
#include <iomanip>
#include <time.h>

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

AbstractOptimizer::AbstractOptimizer(Function * function):
  _function(function) 
{
	// Initialization with defaults:
	_messageHandler   = & cout;
	_profiler         = & cout;
	_constraintPolicy = CONSTRAINTS_KEEP;	
	_nbEvalMax        = 1000000;
  _nbEval           = 0;
	_verbose          = true;
  _isInitialized    = false;
}

/******************************************************************************/

AbstractOptimizer::AbstractOptimizer(const AbstractOptimizer & opt)
{
  _function               = opt._function;
  _parameters             = opt._parameters;
  _messageHandler         = opt._messageHandler;
  _profiler               = opt._profiler;
  _constraintPolicy       = opt._constraintPolicy;
  _tolIsReached           = opt._tolIsReached;
  if(opt._stopCondition)
  {
    _stopCondition        = opt._stopCondition->clone();
    _stopCondition->setOptimizer(this);
  }
  else
    _stopCondition        = NULL;
  if(opt._defaultStopCondition)
  {
    _defaultStopCondition = opt._defaultStopCondition->clone();
    _defaultStopCondition->setOptimizer(this);
  }
  else
    _defaultStopCondition = NULL;
  _nbEvalMax              = opt._nbEvalMax;
  _nbEval                 = opt._nbEval;
  _verbose                = opt._verbose;
  //In case of AutoParameter instances, we must actualize the pointers toward _messageHandler:
  init(_parameters);
  _isInitialized          = opt._isInitialized;
}

/******************************************************************************/

AbstractOptimizer & AbstractOptimizer::operator=(const AbstractOptimizer & opt)
{
  _function               = opt._function;
  _parameters             = opt._parameters;
  _messageHandler         = opt._messageHandler;
  _profiler               = opt._profiler;
  _constraintPolicy       = opt._constraintPolicy;
  _tolIsReached           = opt._tolIsReached;
  if(opt._stopCondition)
  {
    _stopCondition        = opt._stopCondition->clone();
    _stopCondition->setOptimizer(this);
  }
  else
    _stopCondition        = NULL;
  if(opt._defaultStopCondition)
  {
    _defaultStopCondition = opt._defaultStopCondition->clone();
    _defaultStopCondition->setOptimizer(this);
  }
  else
    _defaultStopCondition = NULL;
  _nbEvalMax              = opt._nbEvalMax;
  _nbEval                 = opt._nbEval;
  _verbose                = opt._verbose;
  //In case of AutoParameter instances, we must actualize the pointers toward _messageHandler:
  init(_parameters);
  _isInitialized          = opt._isInitialized;
  return *this;
}

/******************************************************************************/
	
void AbstractOptimizer::init(const ParameterList & params) throw (Exception)
{
	_parameters = params;
	     if(_constraintPolicy == CONSTRAINTS_AUTO)   autoParameter();
	else if(_constraintPolicy == CONSTRAINTS_IGNORE) ignoreConstraints();
	_tolIsReached = false;
  _isInitialized = true;
}

/******************************************************************************/

void AbstractOptimizer::profile(double v)
{
	if(_profiler != NULL) (* _profiler) << v;
}	

/******************************************************************************/

void AbstractOptimizer::profileln(double v)
{
	if(_profiler != NULL) (* _profiler) << v << endl;
}
	
/******************************************************************************/

void AbstractOptimizer::profile(const string & s)
{
	if(_profiler != NULL) (* _profiler) << s;
}	

/******************************************************************************/

void AbstractOptimizer::profileln(const string & s)
{
	if(_profiler != NULL) (* _profiler) << s << endl;
}
	
/******************************************************************************/

void AbstractOptimizer::printPoint(const ParameterList & params, double value)
{
	int ndim = params.size();
	for (int j = 0; j < ndim; j++)
  {
		profile(TextTools::toString(params[j]->getValue()));
		profile("\t"); 
	}
	profile(value);
  profile("\t");
  time_t seconds;
  time(&seconds);
  profileln(ctime(&seconds)); 
}

/******************************************************************************/

void AbstractOptimizer::printMessage(const string & message)
{
	if(_messageHandler != NULL) (* _messageHandler) << message << endl;
}

/******************************************************************************/

void AbstractOptimizer::autoParameter()
{
	for(unsigned int i = 0; i < _parameters.size(); i++)
  {
		Parameter * p = _parameters[i];
		AutoParameter * ap = new AutoParameter(* p);
		ap->setMessageHandler(_messageHandler);
		_parameters[i] = ap;
		delete p;
	}
}

/******************************************************************************/

void AbstractOptimizer::ignoreConstraints()
{
	for(unsigned int i = 0; i < _parameters.size(); i++) {
		_parameters[i] -> removeConstraint();
	}
}

/******************************************************************************/

string AbstractOptimizer::CONSTRAINTS_AUTO   = "auto";
string AbstractOptimizer::CONSTRAINTS_IGNORE = "ignore";
string AbstractOptimizer::CONSTRAINTS_KEEP   = "keep";

/******************************************************************************/

