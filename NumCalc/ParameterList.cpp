//
// File: ParameterList.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Oct 15 18:17:29 2003
//

/*
Copyright ou © ou Copr. Julien Dutheil, (19 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour le calcul numérique.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. Julien Dutheil, (November 19, 2004)

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
		
ParameterList & ParameterList::operator=(const ParameterList & pl)
{
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
	
ParameterList::~ParameterList()
{
	// Delete all parameter objects.
	reset();
}
	
/******************************************************************************/

inline const Parameter * ParameterList::getParameter(const string & name) const
{
	for(unsigned int i = 0; i < size(); i++) {
		const Parameter * p = operator[](i);
		if(p -> getName() == name) return p;
	}
	return NULL;
}

/******************************************************************************/

inline Parameter * ParameterList::getParameter(const string & name)
{
	for(unsigned int i = 0; i < size(); i++) {
		Parameter * p = operator[](i);
		if(p -> getName() == name) return p;
	}
	return NULL;
}

/******************************************************************************/

inline ParameterList ParameterList::subList(const vector<string> & names) const
{
	ParameterList pl;
	for(unsigned int i = 0; i < names.size(); i++) {
		const Parameter * param = getParameter(names[i]);
		if(param != NULL) pl.push_back(dynamic_cast<Parameter *>(param -> clone()));
	}
	return pl;
}

/******************************************************************************/
inline ParameterList ParameterList::subList(const string & name) const
{
	ParameterList pl;
	const Parameter * param = getParameter(name);
	if(param != NULL) pl.push_back(dynamic_cast<Parameter *>(param -> clone()));
	return pl;
}

/******************************************************************************/

inline ParameterList ParameterList::subList(vector<unsigned int> parameters) const 
{
	ParameterList pl;
	for(unsigned int i = 0; i < parameters.size(); i++) {
		if(parameters[i] < size()) pl.push_back(dynamic_cast<Parameter *>(this -> operator[](parameters[i]) -> clone()));
	}
	return pl;
}

/******************************************************************************/

inline ParameterList ParameterList::subList(unsigned int parameter) const
{
	ParameterList pl;
  if(parameter < size()) pl.push_back(dynamic_cast<Parameter *>(this -> operator[](parameter) -> clone()));
	return pl;
}

/******************************************************************************/
	
inline ParameterList ParameterList::getCommonParametersWith(const ParameterList & params) const
{
	ParameterList pl;
  for(unsigned int i = 0; i < params.size(); i++) {
		Parameter * p = params[i];
		if(getParameter(p -> getName()) != NULL)
			pl.push_back(dynamic_cast<Parameter *>(p -> clone()));
	}
	
	return pl;
}
	
/******************************************************************************/

inline vector<string> ParameterList::getParameterNames() const
{
	vector<string> pNames(size());
	for(unsigned int i = 0; i < size(); i++) pNames[i] = operator[](i) -> getName();
	return pNames;
}

/******************************************************************************/

void ParameterList::addParameter(const Parameter & param) throw (ParameterException)
{
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

inline void ParameterList::setParameterValue(const string & name, double value)
throw (ParameterNotFoundException, ConstraintException)
{
	Parameter * p = getParameter(name);
	if(p == NULL) throw ParameterNotFoundException("ParameterList::setParameter", name);
	p -> setValue(value);
}
	
/******************************************************************************/

inline void ParameterList::setAllParametersValues(const ParameterList & params)
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

inline void ParameterList::setParametersValues(const ParameterList & params) 
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

inline void ParameterList::matchParametersValues(const ParameterList & params) 
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

inline void ParameterList::setAllParameters(const ParameterList & params)
throw (ParameterNotFoundException)
{
	for(ParameterList::iterator it = begin(); it < end(); it++) {
		const Parameter * p = params.getParameter((*it) -> getName());
		if(p == NULL) throw ParameterNotFoundException("ParameterList::setAllParameters", (* it) -> getName());
		**it = *p;
	}		
}

/******************************************************************************/

inline void ParameterList::setParameters(const ParameterList & params) 
throw (ParameterNotFoundException)
{
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		if(p == NULL) throw ParameterNotFoundException("ParameterList::setParameters", (* it) -> getName());
		*p = **it;
	}
}

/******************************************************************************/

inline void ParameterList::matchParameters(const ParameterList & params) 
{
	for(ParameterList::const_iterator it = params.begin(); it < params.end(); it++) {
		Parameter * p = getParameter((*it) -> getName());
		if(p != NULL) *p = **it;
	}		
}

/******************************************************************************/

inline void ParameterList::deleteParameter(const string & name) throw (ParameterNotFoundException)
{
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

void ParameterList::printParameters(ostream & out) const
{
	out << "Name:\tValue:\tConstraint:" << endl;
	out << "_________________________________________________" << endl;
	for(unsigned int i = 0; i < size(); i++) {
		out << operator[](i) -> getName();
		out	<< "\t" << operator[](i) -> getValue();
		out	<< (operator[](i) -> hasConstraint() ? "\t" + operator[](i) -> getConstraint() -> getDescription() : string(""));
		out << endl;
	}
}

/******************************************************************************/

void ParameterList::reset()
{
	for(unsigned int i = 0; i < size(); i++) {
		//cout << "Deleting parameter " << i << ": " << operator[](i) -> getName() << endl;
		delete operator[](i);
	}
	resize(0);
}

/******************************************************************************/
