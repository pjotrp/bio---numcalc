//
// File: SimpleMultiDimensions.cpp
// Created by; jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: ue Nov 16 17:51 2004
//

/*
Copyright ou � ou Copr. Julien Dutheil, (19 Novembre 2004) 

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
Copyright or � or Copr. Julien Dutheil, (November 19, 2004)

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

/******************************************************************************/

#include "SimpleMultiDimensions.h"

/******************************************************************************/

SimpleMultiDimensions::SimpleMultiDimensions(Function * function):
	AbstractOptimizer(function)
{
	_defaultStopCondition = new FunctionStopCondition(this);
	_stopCondition = _defaultStopCondition;
	_nbParams = 0;
}

/******************************************************************************/

SimpleMultiDimensions::~SimpleMultiDimensions()
{
	for(unsigned int i = 0; i < _nbParams; i++) {
		delete _optimizers[i];
	}
	delete _defaultStopCondition;
}

/******************************************************************************/

void SimpleMultiDimensions::init(const ParameterList & params) throw (Exception)
{
	// Some cleaning first.
	// This is useful only if the optimizer have been initialized once before this time.
	for(unsigned int i = 0; i < _nbParams; i++) {
		delete _optimizers[i];
	}

	_parameters = params;

	_nbParams = params.size();
	if(_nbParams == 0) return;

	// Initialize optimizers:
	unsigned int nbEvalMax = _nbEvalMax / _nbParams;
	_optimizers.resize(_nbParams);
	for(unsigned int i = 0; i < _nbParams; i++) {
		//_optimizers[i] = new BrentOneDimension(_function);
		_optimizers[i] = new GoldenSectionSearch(_function);
		_optimizers[i] -> setMaximumNumberOfEvaluations(nbEvalMax);
		_optimizers[i] -> setProfiler(_profiler);
		_optimizers[i] -> setMessageHandler(_messageHandler);
		//_optimizers[i] -> setTolerance(_stopCondition -> getTolerance());
		_optimizers[i] -> getStopCondition() -> setTolerance(getStopCondition() -> getTolerance());
		_optimizers[i] -> setConstraintPolicy(_constraintPolicy);
		_optimizers[i] -> setInitialInterval(0.,1.);
		profile(_parameters[i] -> getName() + "\t"); 
	}
	
	profileln("Function");

	printPoint(_parameters, _function -> f(_parameters));
	// Initialize the StopCondition:
	_stopCondition -> isToleranceReached();
}

/******************************************************************************/

double SimpleMultiDimensions::step() throw (Exception)
{
	for(unsigned int i = 0; i < _nbParams; i++) {
		// Re-init optimizer according to new values:
		double v = _parameters[i] -> getValue();
		_optimizers[i] -> setInitialInterval(v - 0.01, v + 0.01);
		_optimizers[i] -> init(_parameters.subList(i));
		cout << _parameters.subList(i)[0] -> getValue() << endl;

		// Optimize through this dimension:
		_optimizers[i] -> optimize();
		// Update parameters with the new value:
		cout << _optimizers[i] -> getParameters()[0] -> getValue() << endl;
		_parameters.setParametersValues(_optimizers[i] -> getParameters());
		cout << _parameters[0] -> getValue() << endl;
		_nbEval += _optimizers[i] -> getNumberOfEvaluations(); 
	}
	_tolIsReached = _nbParams <= 1 || _stopCondition -> isToleranceReached();
}

/******************************************************************************/

double SimpleMultiDimensions::optimize() throw (Exception)
{
	_tolIsReached = false;
	_nbEval = 0;
	for (_nbEval = 0; _nbEval < _nbEvalMax && ! _tolIsReached; _nbEval++) {
		cout << "#" << endl;
		step();
	}
	return _function -> getValue();
}

/******************************************************************************/

double SimpleMultiDimensions::getFunctionValue() const
{
 return _function -> getValue();
}

/******************************************************************************/

