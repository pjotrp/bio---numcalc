//
// File: OptimizationStopCondition.cpp
// Created by: Julien Dutheil
// Created on: Tue Dec 23 11:51:31 2003
//

/*
Copyright ou � ou Copr. CNRS, (19 Novembre 2004) 

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

#include "OptimizationStopCondition.h"
#include "Optimizer.h"
#include "VectorTools.h"
using namespace VectorFunctions;
#include "NumTools.h"
using namespace NumTools;

/******************************************************************************/

OptimizationStopCondition::OptimizationStopCondition() {}

OptimizationStopCondition::~OptimizationStopCondition() {}

/******************************************************************************/

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(const Optimizer * optimizer):
	_optimizer(optimizer),
	_tolerance(0.000001),
	_callCount(0),
	_burnin(0) {}

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(
	const Optimizer * optimizer,
	double tolerance):
	_optimizer(optimizer),
	_tolerance(tolerance),
	_callCount(0),
	_burnin(0) {}

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(
	const Optimizer * optimizer,
	int burnin):
	_optimizer(optimizer),
	_tolerance(0.000001),
	_callCount(0),
	_burnin(burnin) {}

AbstractOptimizationStopCondition::AbstractOptimizationStopCondition(
	const Optimizer * optimizer,
	double tolerance,
	int burnin):
	_optimizer(optimizer),
	_tolerance(tolerance),
	_callCount(0),
	_burnin(burnin) {}

AbstractOptimizationStopCondition::~AbstractOptimizationStopCondition() {}

void AbstractOptimizationStopCondition::setTolerance(double tolerance)
{
	_tolerance = tolerance;
}

double AbstractOptimizationStopCondition::getTolerance() const
{
	return _tolerance;
}

void AbstractOptimizationStopCondition::resetCounter()
{
	_callCount = 0;
}

void AbstractOptimizationStopCondition::setBurnin(int burnin)
{
	_burnin = burnin;
}

int AbstractOptimizationStopCondition::getBurnin() const
{
	return _burnin;
}

/******************************************************************************/

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer):
	AbstractOptimizationStopCondition(optimizer)
{
	_newParametersEstimates = _optimizer -> getParameters();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer,
	double tolerance):
	AbstractOptimizationStopCondition(optimizer, tolerance)
{
	_newParametersEstimates = _optimizer -> getParameters();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer,
	int burnin):
	AbstractOptimizationStopCondition(optimizer, burnin)
{
	init();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::ParametersStopCondition(
	const Optimizer * optimizer,
	double tolerance,
	int burnin):
	AbstractOptimizationStopCondition(optimizer, tolerance, burnin)
{
	init();
	if(_newParametersEstimates.size() == 0) {
		cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
		     << "Be sure to have initialized the Optimizer first!" << endl; 
	}
}

ParametersStopCondition::~ParametersStopCondition() {}

/******************************************************************************/

void ParametersStopCondition::init()
{
	if(_optimizer -> getFunction() != NULL)
		_newParametersEstimates = _optimizer -> getParameters();
}

/******************************************************************************/

bool ParametersStopCondition::isToleranceReached() const
{
	_callCount++;
	_lastParametersEstimates = _newParametersEstimates;
	_newParametersEstimates   = _optimizer -> getParameters();
	if(_callCount <= _burnin) return false;
	for(unsigned int i = 0; i < _newParametersEstimates.size(); i++) {
		Parameter * p = _newParametersEstimates[i];
		if(p == NULL) throw ParameterNotFoundException("ParameterStopCondition::isToleranceReached.", p -> getName());
		double lastEstimate = _lastParametersEstimates.getParameter(p -> getName()) -> getValue();
		double newEstimate = p -> getValue();
		double tol = NumTools::abs<double>(newEstimate - lastEstimate);
		if(tol > _tolerance) {
			return false;
		}
	}
	return true;
}

/******************************************************************************/

FunctionStopCondition::FunctionStopCondition(
	const Optimizer * optimizer):
	AbstractOptimizationStopCondition(optimizer)
{
	init();
}

FunctionStopCondition::FunctionStopCondition(
	const Optimizer * optimizer,
	double tolerance):
	AbstractOptimizationStopCondition(optimizer, tolerance)
{
	init();
}

FunctionStopCondition::FunctionStopCondition(
	const Optimizer * optimizer,
	int burnin):
	AbstractOptimizationStopCondition(optimizer, burnin)
{
	init();
}

FunctionStopCondition::FunctionStopCondition(
	const Optimizer * optimizer,
	double tolerance,
	int burnin):
	AbstractOptimizationStopCondition(optimizer, tolerance, burnin)
{
	init();
}

FunctionStopCondition::~FunctionStopCondition() {}

/******************************************************************************/

void FunctionStopCondition::init()
{
	if(_optimizer -> getFunction() != NULL)
		_newFunctionValue = _optimizer -> getFunctionValue();
}

/******************************************************************************/

bool FunctionStopCondition::isToleranceReached() const
{
	_callCount++;
	_lastFunctionValue = _newFunctionValue;
	_newFunctionValue  = _optimizer -> getFunctionValue();
	if(_callCount <= _burnin) return false;
	double tol = NumTools::abs<double>(_newFunctionValue - _lastFunctionValue);
	return tol < _tolerance;
}

/******************************************************************************/
