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

SimpleMultiDimensions::SimpleMultiDimensions(const Function * function):
	AbstractOptimizer(function)
{
	_stopCondition = new FunctionStopCondition(this); 
}

/******************************************************************************/

SimpleMultiDimensions::~SimpleMultiDimensions()
{
	for(unsigned int i = 0; i < _nbParams; i++) {
		delete _optimizers[i];
	}
	delete _stopCondition;
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
	// Initialize optimizers:
	unsigned int nbEvalMax = (unsigned int)(_nbEvalMax / _nbParams);
	_optimizers.resize(_nbParams);
	for(unsigned int i = 0; i < _nbParams; i++) {
		_optimizers[i] = new BrentOneDimension(_function);
		_optimizers[i] -> setMaximumNumberOfEvaluations(nbEvalMax);
	}	
}

/******************************************************************************/

double SimpleMultiDimensions::step() throw (Exception)
{
	for(unsigned int i = 0; i < _nbParams; i++) {
		// Re-init optimizer according to new values:
		_optimizers[i] -> init(_parameters.subList(i));
		// Optimize through this dimension:
		_optimizers[i] -> optimize();
		// Update parameters with the new value:
		_parameters.setParametersValues(_optimizers[i] -> getParameters());
	}
	_tolIsReached = _stopCondition -> isToleranceReached();
}

/******************************************************************************/

double SimpleMultiDimensions::optimize() throw (Exception)
{
	_tolIsReached = false;
	for (_nbEval = 0; _nbEval < _nbEvalMax && ! _tolIsReached; _nbEval++) {
		step();
	}
	return _currentValue;
}

/******************************************************************************/

