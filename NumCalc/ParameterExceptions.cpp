//
// File: ParameterExceptions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov  3 18:05:36 2003
//

#include "ParameterExceptions.h"
#include "Parameter.h"

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

ParameterException::ParameterException(const char *   text, const Parameter * param) :
	Exception("ParameterException: " + string(text) + (param != NULL ? "(" + param -> getName() + ")" : string(""))),
	parameter(param) {};
		
ParameterException::ParameterException(const string & text, const Parameter * param) :
	Exception("ParameterException: " + text + (param != NULL ? "(" + param -> getName() + ")" : string(""))),
	parameter(param) {};
		
ParameterException::~ParameterException() throw() {};
	
const Parameter * ParameterException::getParameter() const { return parameter; }

/******************************************************************************/

ConstraintException::ConstraintException(const char *   text, const Parameter * param, double badValue) :
	ParameterException("ConstraintException: " + string(text) + "(" + TextTools::toString(badValue) + ")", param),
	badValue(badValue) {};
		
ConstraintException::ConstraintException(const string & text, const Parameter * param, double badValue) :
	ParameterException("ConstraintException: " + text + "(" + TextTools::toString(badValue) + ")", param),
	badValue(badValue) {};
		
ConstraintException::~ConstraintException() throw() {};
	
double ConstraintException::getBadValue() const { return badValue; }

/******************************************************************************/

ParameterNotFoundException::ParameterNotFoundException(const char *   text, const string & param) :
	Exception("ParameterNotFoundException: " + string(text) + (!TextTools::isEmpty(param) ? "(" + param + ")" : string(""))),
	parameter(param) {};
		
ParameterNotFoundException::ParameterNotFoundException(const string & text, const string & param) :
	Exception("ParameternotFoundException: " + text + (!TextTools::isEmpty(param) ? "(" + param + ")" : string(""))),
	parameter(param) {};
		
ParameterNotFoundException::~ParameterNotFoundException() throw() {};
	
string ParameterNotFoundException::getParameter() const { return parameter; }

/******************************************************************************/
