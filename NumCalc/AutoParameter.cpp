//
// File: AutoParameter.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov 11 22:15:16 2003
//

#include "AutoParameter.h"

#include <iostream>

using namespace std;

// Utils:
#include "Utils/TextTools.h"

double AutoParameter::TINY = .000000000001;

/** Constructors: *************************************************************/

AutoParameter::AutoParameter(const string & name, double value, const Constraint * constraint):
Parameter(name, value, constraint)
{
	_messageHandler = &cout;
}

AutoParameter::AutoParameter(const Parameter & p): Parameter(p)
{
	_messageHandler = &cout;
}

AutoParameter::AutoParameter(const AutoParameter & p): Parameter(p)
{
	_messageHandler = p._messageHandler;
}

AutoParameter & AutoParameter::operator=(const Parameter & p) {
	name       = p.getName();
	value      = p.getValue();
	constraint = p.getConstraint();
	try {
		_messageHandler = dynamic_cast<const AutoParameter &>(p)._messageHandler;
	} catch(exception e) {
		_messageHandler = &cout;
	}
	return * this;	
}

/** Destructor: ***************************************************************/

AutoParameter::~AutoParameter() {}

/** Clonable interface: *******************************************************/

Clonable * AutoParameter::clone() const { return new AutoParameter(* this); }

/******************************************************************************/
	
void AutoParameter::setValue(double value) throw (ConstraintException) {
	try { // First we try to assign this value:
		Parameter::setValue(value);
	} catch (ConstraintException ce) { // Aie, there's a pb here...
		if(_messageHandler != NULL) {
			(* _messageHandler) << "Constraint match at parameter ";
			(* _messageHandler) << name;
			(* _messageHandler) << ", badValue = ";
			(* _messageHandler) << ce.getBadValue();
			(* _messageHandler) << endl;
		}
		double limit = constraint -> getLimit(value);
		try { // We try to assign the limit then.
			Parameter::setValue(limit);
		} catch(ConstraintException ce2) { // Aie, the limit is not reachable, so we perform a smaller step...
			//Parameter::setValue((getValue() + limit) / 2);
			try {
				// Try on the right:
				Parameter::setValue(limit + TINY);
			} catch(ConstraintException ce3) {
				// Try on the left:
				Parameter::setValue(limit - TINY);
			}
		}
	}
}

/******************************************************************************/

void AutoParameter::setMessageHandler(ostream * mh) { _messageHandler = mh; }

/******************************************************************************/
