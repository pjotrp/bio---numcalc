//
// File: ParameterList.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Oct 15 18:17:29 2003
//

#include "ParameterList.h"

#include <iostream>

/** Constructors: *************************************************************/

ParameterList::ParameterList() : vector<Parameter *>() {}

/** Copy constructor: *********************************************************/
	
ParameterList::ParameterList(const ParameterList & pl): vector<Parameter*>(pl.size()) {
	// First delete all parameters:
	//reset();
	
	// Then resize the vector:
	//resize(pl.size());
	
	// Now copy all parameters:
	for(unsigned int i = 0; i < size(); i++) {
		operator[](i) = dynamic_cast<Parameter *>(pl[i] -> clone());
	}
}

/** Assignation operator: *****************************************************/
		
ParameterList & ParameterList::operator=(const ParameterList & pl) {
	// First delete all parameters:
	reset();
		
	// Then resize the vector:
	resize(pl.size());
	
	// Now copy all parameters:
	for(unsigned int i = 0; i < size(); i++) {
		operator[](i) = dynamic_cast<Parameter *>(pl[i] -> clone());
	}
	
	return * this;
}
	
/** Destructor: ***************************************************************/
	
ParameterList::~ParameterList() {
	// Delete all parameter objects.
	reset();
}
	
/******************************************************************************/

const Parameter * ParameterList::getParameter(const string & name) const {
	for(unsigned int i = 0; i < size(); i++) {
		const Parameter * p = operator[](i);
		if(p -> getName() == name) return p;
	}
	return NULL;
}

/******************************************************************************/

Parameter * ParameterList::getParameter(const string & name) {
	for(unsigned int i = 0; i < size(); i++) {
		Parameter * p = operator[](i);
		if(p -> getName() == name) return p;
	}
	return NULL;
}

/******************************************************************************/

vector<string> ParameterList::getParameterNames() const {
	vector<string> pNames(size());
	for(unsigned int i = 0; i < size(); i++) pNames[i] = operator[](i) -> getName();
	return pNames;
}

/******************************************************************************/

void ParameterList::addParameter(const Parameter & param) throw (ParameterException) {
	if(getParameter(param.getName()) != NULL)
		throw ParameterException("ParameterList::addParameter. Parameter with name '" + param.getName() + "' already exists.", &param); 
	push_back(dynamic_cast<Parameter *>(param.clone()));
}

/******************************************************************************/

void ParameterList::addParameters(const ParameterList & params)
throw (ParameterException)
{
	for(unsigned int i = 0; i < params.size(); i++) addParameter(* params[i]);
}

/******************************************************************************/

void ParameterList::setParameterValue(const string & name, double value)
throw (ParameterNotFoundException, ConstraintException)
{
	Parameter * p = getParameter(name);
	if(p == NULL) throw ParameterNotFoundException("ParameterList::setParameter", name);
	p -> setValue(value);
}
	
/******************************************************************************/

void ParameterList::setAllParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{
	// First we check if all values are correct:
	for(ParameterList::iterator it = begin(); it < end(); it++) {
		const Parameter * p = params.getParameter((*it) -> getName());
		if(p == NULL) throw ParameterNotFoundException("ParameterList::setAllParameters", (* it) -> getName());
		if((*it) -> hasConstraint() && !(*it) -> getConstraint() -> isCorrect(p -> getValue()))
			throw ConstraintException("ParameterList::setParametersValues()", *it, p -> getValue());
	}		

	// If all values are ok, we set them: 
	for(ParameterList::iterator it = begin(); it < end(); it++) {
		const Parameter * p = params.getParameter((*it) -> getName());
		(*it) -> setValue(p -> getValue());
	}		
}

/******************************************************************************/

void ParameterList::setParametersValues(const ParameterList & params) 
throw (ParameterNotFoundException, ConstraintException)
{
	// First we check if all values are correct:
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		if(p == NULL) throw ParameterNotFoundException("ParameterList::setParameters", (* it) -> getName());
		if(p -> hasConstraint() && !p -> getConstraint() -> isCorrect((*it) -> getValue()))
			throw ConstraintException("ParameterList::setParametersValues()", p, (*it) -> getValue());
	}
	
	// If all values are ok, we set them: 
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		p -> setValue((*it) -> getValue());
	}
}

/******************************************************************************/

void ParameterList::matchParametersValues(const ParameterList & params) 
throw (ConstraintException) 
{
	// First we check if all values are correct:
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		if(p != NULL) if(p -> hasConstraint() && !p -> getConstraint() -> isCorrect((*it) -> getValue()))
			throw ConstraintException("ParameterList::mathcParametersValues()", p, (*it) -> getValue());
	}		

	// If all values are ok, we set them: 
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		if(p != NULL) p -> setValue((*it) -> getValue());
	}		
}

/******************************************************************************/

void ParameterList::setAllParameters(const ParameterList & params)
throw (ParameterNotFoundException)
{
	for(ParameterList::iterator it = begin(); it < end(); it++) {
		const Parameter * p = params.getParameter((*it) -> getName());
		if(p == NULL) throw ParameterNotFoundException("ParameterList::setAllParameters", (* it) -> getName());
		**it = *p;
	}		
}

/******************************************************************************/

void ParameterList::setParameters(const ParameterList & params) 
throw (ParameterNotFoundException)
{
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		if(p == NULL) throw ParameterNotFoundException("ParameterList::setParameters", (* it) -> getName());
		*p = **it;
	}
}

/******************************************************************************/

void ParameterList::matchParameters(const ParameterList & params) 
{
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		if(p != NULL) *p = **it;
	}		
}

/******************************************************************************/

void ParameterList::deleteParameter(const string & name) throw (ParameterNotFoundException) {
	for(unsigned int i = 0; i < size(); i++) {
		Parameter * p = operator[](i);
		if(p -> getName() == name) {
			delete p;
			erase(begin() + i);
			return;
		}
	}
	throw ParameterNotFoundException("ParameterList::deleteParameter", name);
}

/******************************************************************************/

void ParameterList::printParameters(ostream & out) const {
	out << "Name:\tValue:\tConstraint:" << endl;
	out << "_________________________________________________" << endl;
	for(unsigned int i = 0; i < size(); i++) {
		out << operator[](i) -> getName();
		out	<< "\t" << operator[](i) -> getValue();
		out	<< (operator[](i) -> hasConstraint() ? "\t" + operator[](i) -> getConstraint() -> getDescription() : "");
		out << endl;
	}
}

/******************************************************************************/

void ParameterList::reset() {
	for(unsigned int i = 0; i < size(); i++) {
		//cout << "Deleting parameter " << i << ": " << operator[](i) -> getName() << endl;
		delete operator[](i);
	}
}

/******************************************************************************/
