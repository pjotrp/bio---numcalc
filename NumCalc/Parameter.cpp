//
// File: Parameter.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Oct 15 15:40:47 2003
//

#include "Parameter.h"

//From Utils:
#include <Utils/TextTools.h>

#include <iostream>

/** Constructors: *************************************************************/

Parameter::Parameter(const string & name, double value, const Constraint * constraint) 
throw (ConstraintException) {
	this -> name = name;
	this -> constraint = constraint;
	// This may throw a ConstraintException:
	setValue(value);
}

Parameter::Parameter(const Parameter & p): Clonable() {
	name       = p.getName();
	value      = p.getValue();
	constraint = p.getConstraint();	
}

Parameter & Parameter::operator=(const Parameter & p) {
	name       = p.getName();
	value      = p.getValue();
	constraint = p.getConstraint();
	return * this;	
}

/** Destructor: ***************************************************************/

Parameter::~Parameter() {}

/** Clonable interface: *******************************************************/

Clonable * Parameter::clone() const { return new Parameter(* this); }
	
/** Name: *********************************************************************/
	
void Parameter::setName(const string & name) { this -> name = name;}

string Parameter::getName() const {	return name; }

/** Value: ********************************************************************/

void Parameter::setValue(double value) throw (ConstraintException) {
	if(constraint != NULL && !constraint -> isCorrect(value)) 
		throw ConstraintException("Parameter::setValue", this, value);
	this -> value = value;
}

double Parameter::getValue() const { return value; }

/** Constraint: ***************************************************************/

const Constraint * Parameter::getConstraint() const { return constraint; }

bool Parameter::hasConstraint() const { return (constraint != NULL); }

const Constraint * Parameter::removeConstraint() {
	const Constraint * c = constraint;
	constraint = NULL;
	return c;
}

/******************************************************************************/

IncludingPositiveReal Parameter::R_PLUS(0);

ExcludingPositiveReal Parameter::R_PLUS_STAR(0);

/******************************************************************************/
