//
// File: AbstractOptimizer.cpp
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Dec 22 12:18:09 2003
//

/*
Copyright ou � ou Copr. CNRS, (17 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant � fournir des classes
pour le calcul num�rique.

Ce logiciel est r�gi par la licence CeCILL soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffus�e par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant 
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe � 
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement, 
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�. 

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accept� les
termes.
*/

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

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

#include "AbstractOptimizer.h"
#include "AutoParameter.h"

// From the STL:
#include <iomanip>

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

AbstractOptimizer::AbstractOptimizer(Function * function): _function(function) 
{
	// Initialization with defaults:
	_messageHandler   = & cout;
	_profiler         = & cout;
	_constraintPolicy = CONSTRAINTS_KEEP;	
	_nbEvalMax = 1000000;
	_verbose = true;
}

/******************************************************************************/

AbstractOptimizer::~AbstractOptimizer() {}
	
/******************************************************************************/
	
void AbstractOptimizer::init(const ParameterList & params) throw (Exception) {
	_parameters = params;
	     if(_constraintPolicy == CONSTRAINTS_AUTO)   autoParameter();
	else if(_constraintPolicy == CONSTRAINTS_IGNORE) ignoreConstraints();
	_tolIsReached = false;
}

/******************************************************************************/

ParameterList AbstractOptimizer::getParameters() const { return _parameters; }

/******************************************************************************/

void AbstractOptimizer::setMessageHandler(ostream * mh) { _messageHandler = mh; }

/******************************************************************************/

void AbstractOptimizer::setProfiler(ostream * profiler) { _profiler = profiler; }

/******************************************************************************/

void AbstractOptimizer::setStopCondition(
	OptimizationStopCondition * stopCondition)
{
	_stopCondition = stopCondition;
}

/******************************************************************************/

OptimizationStopCondition * AbstractOptimizer::getStopCondition()
{
	return _stopCondition;
}

/******************************************************************************/

OptimizationStopCondition * AbstractOptimizer::getDefaultStopCondition()
{
	return _defaultStopCondition;
}

/******************************************************************************/

bool AbstractOptimizer::isToleranceReached() const
{
	return _tolIsReached;
}

/******************************************************************************/

bool AbstractOptimizer::isMaximumNumberOfEvaluationsReached() const {
	return _nbEval >= _nbEvalMax;
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
	for (int j = 0; j < ndim; j++) {
		profile(TextTools::toString(params[j] -> getValue()));
		profile("\t"); 
	}
	profileln(value);
}

/******************************************************************************/

void AbstractOptimizer::printMessage(const string & message)
{
	if(_messageHandler != NULL) (* _messageHandler) << message << endl;
}

/******************************************************************************/

void AbstractOptimizer::autoParameter()
{
	for(unsigned int i = 0; i < _parameters.size(); i++) {
		Parameter * p = _parameters[i];
		AutoParameter * ap = new AutoParameter(* p);
		ap -> setMessageHandler(_messageHandler);
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

/** Constraint policy: ********************************************************/

void AbstractOptimizer::setConstraintPolicy(const string & constraintPolicy)
{ _constraintPolicy = constraintPolicy; }

string AbstractOptimizer::getConstraintPolicy() const { return _constraintPolicy; }

string AbstractOptimizer::CONSTRAINTS_AUTO   = "auto";
string AbstractOptimizer::CONSTRAINTS_IGNORE = "ignore";
string AbstractOptimizer::CONSTRAINTS_KEEP   = "keep";

/******************************************************************************/

void AbstractOptimizer::setMaximumNumberOfEvaluations(int max) { _nbEvalMax = max; }

/******************************************************************************/

int AbstractOptimizer::getNumberOfEvaluations() const { return _nbEval; }

/******************************************************************************/

