//
// File: BrentOneDimension.cpp
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 17 11:45:58 2003
//

/*
Copyright ou © ou Copr. CNRS, (17 Novembre 2004) 

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
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "BrentOneDimension.h"
#include "NumTools.h"
#include "OneDimensionOptimizationTools.h"

using namespace NumTools;

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

BrentOneDimension::BODStopCondition::BODStopCondition(BrentOneDimension * bod):
	AbstractOptimizationStopCondition(bod) {}

BrentOneDimension::BODStopCondition::~BODStopCondition() {}
		
bool BrentOneDimension::BODStopCondition::isToleranceReached() const
{
	// NRC Test for done:
	const BrentOneDimension * bod = dynamic_cast<const BrentOneDimension *>(_optimizer);
	double x    = bod -> x;
	double xm   = bod -> xm;
	double tol2 = bod -> tol2;
	double b    = bod -> b;
	double a    = bod -> a;
	return NumTools::abs(x - xm) <= (tol2 - 0.5 * (b - a));
}
		
/******************************************************************************/

BrentOneDimension::BrentOneDimension(Function * function) :
AbstractOptimizer(function)
{
	_defaultStopCondition = new BODStopCondition(this);
	_stopCondition = _defaultStopCondition;
	_nbEvalMax = 10000;
}


BrentOneDimension::~BrentOneDimension()
{
	delete _defaultStopCondition;
}
	
/******************************************************************************/
	
double BrentOneDimension::CGOLD = 0.3819660;
double BrentOneDimension::ZEPS  = 1.0e-10;
	
/******************************************************************************/
	
void BrentOneDimension::init(const ParameterList & params) throw (Exception)
{
	// Set the initial value (no use here! Use setInitialValues() instead).
	if(params.size() != 1) throw Exception("BrentOneDimension::init(). This optimizer only deals with one parameter.");
	AbstractOptimizer::init(params);

	// Bracket the minimum.
	Bracket bracket = OneDimensionOptimizationTools::bracketMinimum(_xinf, _xsup, _function, _parameters);
	if(_verbose > 0) {
		printMessage("Initial bracketing:");
		printMessage("A: x = " + TextTools::toString(bracket.a.x) + ", f = " + TextTools::toString(bracket.a.f));
		printMessage("B: x = " + TextTools::toString(bracket.b.x) + ", f = " + TextTools::toString(bracket.b.f));
		printMessage("C: x = " + TextTools::toString(bracket.c.x) + ", f = " + TextTools::toString(bracket.c.f));
	}
	
	// This will be the distance moved on the step before last.
	e = 0.0;

	// a and b must be in ascending order, but input abscissa need not be.
	a = (bracket.a.x < bracket.c.x ? bracket.a.x : bracket.c.x);
	b = (bracket.a.x > bracket.c.x ? bracket.a.x : bracket.c.x);
	
	// Initializations...
	x = w = v = bracket.b.x;
	_parameters[0] -> setValue(x);
	fw = fv = fx = _function -> f(_parameters);	

	profileln(_parameters[0] -> getName() + "\tFunction");
}

/******************************************************************************/
	
void BrentOneDimension::setInitialInterval(double inf, double sup)
{
	_xinf = inf; _xsup = sup;	
}

/******************************************************************************/

void BrentOneDimension::setTolerance(double tolerance)
{
	_tolerance = tolerance;
}

/******************************************************************************/

double BrentOneDimension::step() throw (Exception) {
	if(_verbose > 0) { cout << "."; cout.flush(); }
	xm   = 0.5 * (a + b);
	tol2 = 2.0 * (tol1 = _tolerance * NumTools::abs(x) + ZEPS);
	
	// Construct a trial parabolic fit:
	if (NumTools::abs(e) > tol1) {
		r = (x - w) * (fx - fv);
		q = (x - v) * (fx - fw);
		p = (x - v) * q - (x - w) * r;
		q = 2.0 * (q - r);
		if (q > 0.0) p = -p;
		q = NumTools::abs(q);
		etemp = e;
		e = d;
		if (NumTools::abs(p) >= NumTools::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		// The above conditions determine the acceptability of the parabolic fit.
		// Here we take the golden section step into the larger of the two segments.
		else {
			d = p / q; // Take the parabolic step.
			u = x + d;
			if (u - a < tol2 || b - u < tol2)
				d = sign(tol1, xm - x);
		}
	} else {
		d = CGOLD * (e = (x >= xm ? a - x : b - x));
	}
	u = (NumTools::abs(d) >= tol1 ? x + d : x + sign(tol1, d));

	// This is the one function evaluate per iteration.
	ParameterList pl = _parameters;
	pl[0] -> setValue(u);
	fu = _function -> f(pl);
	// Now decide what to do with our function evaluation.
	if (fu <= fx) {
		if (u >= x) a = x; else b = x;
		// Here is the house keeping:
		shift(v, w, x, u);
		shift(fv, fw, fx, fu);
	} else {
		if (u < x) a = u; else b = u;
		if (fu <= fw || w == x) {
			v  = w;
			w  = u;
			fv = fw;
			fw = fu;
		} else if (fu <= fv || v == x || v == w) {
			v  = u;
			fv = fu;
		}
	} // Done with housekeeping.

	// Store results for this step:
	_parameters[0] -> setValue(x);
	_tolIsReached = _nbEval > 2 && _stopCondition -> isToleranceReached();
	printPoint(_parameters, fx);
	return fx;
}

/******************************************************************************/
	
double BrentOneDimension::optimize() throw (Exception)
{
	// Main program loop.
	_tolIsReached = false;
	for (_nbEval = 0; _nbEval < _nbEvalMax && !_tolIsReached; _nbEval++) {
		step();
	}

	// Apply parameters and evaluate function at minimum point:
	_parameters[0] -> setValue(x); fx = _function -> f(_parameters);
	return fx;
}

/******************************************************************************/

double BrentOneDimension::getFunctionValue() const throw (Exception)
{
	return fx;
}

/******************************************************************************/
