//
// File: OneDimensionOptimizationTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 17 11:15:22 2003
//

#include "OneDimensionOptimizationTools.h"
#include "NumTools.h"

/******************************************************************************
 *                              The Point class                               *
 ******************************************************************************/
 
Point::Point() {}

Point::Point(double x, double f): x(x), f(f) {}

Point::~Point() {}
	
inline void Point::set(double x, double f) { this -> x = x; this -> f = f; }	

/******************************************************************************
 *                             The Bracket class                              *
 ******************************************************************************/
 
Bracket::Bracket() {}

Bracket::~Bracket() {}

/******************************************************************************/

inline void Bracket::setA(double xa, double fa) { a.set(xa, fa); }
inline void Bracket::setB(double xb, double fb) { b.set(xb, fb); }
inline void Bracket::setC(double xc, double fc) { c.set(xc, fc); }

/******************************************************************************/
	
Bracket OneDimensionOptimizationTools::bracketMinimum(
	double a,
	double b,
	Function * function,
	ParameterList parameters)
{
	Bracket bracket;
	// Copy the parameter to use.
	bracket.a.x = a;
	parameters[0] -> setValue(bracket.a.x); bracket.a.f = function -> f(parameters);
	bracket.b.x = b;
	parameters[0] -> setValue(bracket.b.x); bracket.b.f = function -> f(parameters);
	if (bracket.b.f > bracket.a.f) {		
		// Switch roles of first and second point so that we can go downhill
		// in the direction from a to b.
		NumTools::swap<double>(bracket.a.x, bracket.b.x);
		NumTools::swap<double>(bracket.a.f, bracket.b.f);
	}
	
	// First guess for third point:
	bracket.c.x = bracket.b.x + GOLD * (bracket.b.x - bracket.a.x);
	parameters[0] -> setValue(bracket.c.x); bracket.c.f = function -> f(parameters);
	
	// Keep returning here until we bracket:
	while (bracket.b.f > bracket.c.f) {
		// Compute xu by parabolic extrapolation from a, b, c. TINY is used to prevent
		// any possible division by 0.
		double r = (bracket.b.x - bracket.a.x) * (bracket.b.f - bracket.c.f);
		double q = (bracket.b.x - bracket.c.x) * (bracket.b.f - bracket.a.f);
		
		double xu = bracket.b.x - ((bracket.b.x - bracket.c.x) * q - (bracket.b.x - bracket.a.x) * r) /
			(2.0 * NumTools::sign(NumTools::max(NumTools::abs(q - r), TINY), q - r));
		double xulim = (bracket.b.x) + GLIMIT * (bracket.c.x - bracket.b.x);
		double fu;
		
		// We don't go farther than this.
		// Test various possibilities:
		if ((bracket.b.x - xu) * (xu - bracket.c.x) > 0.0) {
			parameters[0] -> setValue(xu); fu = function -> f(parameters);
			if (fu < bracket.c.f) {
				bracket.setA(bracket.b.x, bracket.b.f);
				bracket.setB(xu, fu);
				return bracket;
			} else if (fu > bracket.b.f) {
				bracket.setC(xu, fu);
				return bracket;
			}
			// Parabolic fit was no use.
			// Use default magnification.
			xu = bracket.c.x + GOLD * (bracket.c.x - bracket.b.x);
			parameters[0] -> setValue(xu); fu = function -> f(parameters);
		} else if ((bracket.c.x - xu) * (xu - xulim) > 0.0) {
			// Parabolic fit is between point 3 and its allowed limit.
			parameters[0] -> setValue(xu); fu = function -> f(parameters);
			if (fu < bracket.c.f) {
				NumTools::shift<double>(bracket.b.x, bracket.c.x, xu, bracket.c.x + GOLD * (bracket.c.x - bracket.b.x));
				parameters[0] -> setValue(xu);
				NumTools::shift<double>(bracket.b.f, bracket.c.f, fu, function -> f(parameters));
			}
		} else if ((xu - xulim) * (xulim - bracket.c.x) >= 0.0) {
			// Limit parabolic xu to maximum allowed value.
			xu = xulim;
			parameters[0] -> setValue(xu); fu = function -> f(parameters);
		} else {
			// Reject parabolic xu, use default magnification.
			xu = bracket.c.x + GOLD * (bracket.c.x - bracket.b.x);
			parameters[0] -> setValue(xu); fu = function -> f(parameters);
		}
		// Eliminate oldest point and continue.
		NumTools::shift<double>(bracket.a.x, bracket.b.x, bracket.c.x, xu);
		NumTools::shift<double>(bracket.a.f, bracket.b.f, bracket.c.f, fu);
	}
	return bracket;
}

/******************************************************************************/

double OneDimensionOptimizationTools::GLIMIT = 100.0;
double OneDimensionOptimizationTools::GOLD   = 1.61803399;
double OneDimensionOptimizationTools::TINY   = 1.0e-20;

/******************************************************************************/
