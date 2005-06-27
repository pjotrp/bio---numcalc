//
// File: AbstractDiscreteDistribution.cpp
// Created by: Julien Dutheil
// Created on: ?
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

#include "AbstractDiscreteDistribution.h"

#include "VectorTools.h"
#include "RandomTools.h"

using namespace VectorFunctions;

/******************************************************************************/
	
unsigned int AbstractDiscreteDistribution::getNumberOfCategories() const {
	return _distribution.size();
}

/******************************************************************************/

double AbstractDiscreteDistribution::getCategory(unsigned int categoryIndex) const
{ 
	map<double, double>::const_iterator it = _distribution.begin();
	for(unsigned int i = 0; i < categoryIndex; i++) it++;
	return it -> first;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getProbability(unsigned int categoryIndex) const
{
	map<double, double>::const_iterator it = _distribution.begin();
	for(unsigned int i = 0; i < categoryIndex; i++) it++;
	return it -> second;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getProbability(double category) const 
{
	return _distribution.find(category) -> second;
}

/******************************************************************************/

Vdouble AbstractDiscreteDistribution::getCategories() const
{
	Vdouble result(_distribution.size());
	unsigned int i = 0;
	for(map<double, double>::const_iterator it = _distribution.begin();
		it != _distribution.end();
		it++)
	{
		result[i] = it -> first;
		i++;
	}
	return result;
}

/******************************************************************************/

Vdouble AbstractDiscreteDistribution::getProbabilities() const
{
	Vdouble result(_distribution.size());
	int i = 0;
	for(map<double, double>::const_iterator it = _distribution.begin();
		it != _distribution.end();
		it++) 
	{
		result[i] = it -> second;
		i++;
	}
	return result;
}

/******************************************************************************/

void AbstractDiscreteDistribution::set(double category, double probability) {
	_distribution[category] = probability;
}

/******************************************************************************/

void AbstractDiscreteDistribution::add(double category, double probability) {
	if(_distribution.find(category) == _distribution.end()) {
		//new category
		_distribution[category] = probability;
	} else {
		//existing category
		_distribution[category] += probability;
	}
}

/******************************************************************************/

double AbstractDiscreteDistribution::rand() const 
{
	double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
	double cumprob = 0;
	for(map<double,double>::const_iterator i = _distribution.begin(); 
		i != _distribution.end();
		i++)
	{
		cumprob += i -> second;
		if(r <= cumprob) return i -> first;
	}
	// This line can't be reached:
	return -1.;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getInfCumulativeProbability(double category) const 
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	for(map<double, double>::const_iterator i = _distribution.begin();
		i != it;
		i++) prob += i -> second;
	return prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getIInfCumulativeProbability(double category) const
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	if(it == _distribution.end()) return 0;
	for(map<double, double>::const_iterator i = ++it;
		i != _distribution.end();
		i++) prob += i -> second;
	return 1. - prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getSupCumulativeProbability(double category) const 
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	if(it == _distribution.end()) return 0;
	for(map<double, double>::const_iterator i = ++it;
		i != _distribution.end();
		i++) prob += i -> second;
	return prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getSSupCumulativeProbability(double category) const
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	for(map<double, double>::const_iterator i = _distribution.begin(); 
		i != it;
		i++) prob += i -> second;
	return 1. - prob;
}

/******************************************************************************/

void AbstractDiscreteDistribution::print(ostream & out) const
{
	for(map<double, double>::const_iterator i = _distribution.begin(); i != _distribution.end(); i++) {
		out << "Pr(" << (i -> first) << ") = " << (i -> second) << endl;
	}
}

/******************************************************************************/

