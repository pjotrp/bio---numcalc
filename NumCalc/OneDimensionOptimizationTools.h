//
// File: OneDimensionOptimizationTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 17 11:15:22 2003
//

#ifndef _ONEDIMENSIONOPTIMIZATIONTOOLS_H_
#define _ONEDIMENSIONOPTIMIZATIONTOOLS_H_

#include "Functions.h"

class Point {
			
	public: // Constructor and destructor:
		Point();
		Point(double x, double f);
		virtual ~Point();
		
	public:
		void set(double x, double f);
		
	public:
		double x;
		double f;			
};
		
class Bracket {
	
	public: // Constructor and destructor::
		Bracket();
		virtual ~Bracket();
			
	public: // Methods:
		void setA(double xa, double fa);
		void setB(double xb, double fb);
		void setC(double xc, double fc);
			
	public:
		Point a, b, c;
};


class OneDimensionOptimizationTools
{
	public:
		OneDimensionOptimizationTools() {}
		 ~OneDimensionOptimizationTools() {}

	public:
		
		/**
		 * @brief Bracket a minimum.
		 *
		 * Given a function func, and given distinct initial points x1 and x2,
		 * this routine searches in the downhill direction (defined by the function as
		 * evaluated at the initial points) and returns a Bracket object with new points
		 * a.x, b.x and c.x that bracket a minimum of the function. Also returned are the
		 * function values at the three points, a.f, b.f and c.f.
		 *
		 * @param a, b Two initial values for the parameter.
		 * @return     A bracket object.
		 */
		static Bracket bracketMinimum(double a, double a, Function * function, ParameterList parameters);
	
	public:
		
		/**
		 * @brief Maximum magnification allowed for a parabolic-fit step.
		 */
		static double GLIMIT;
		static double GOLD;
		static double TINY;
	
};


#endif	//_ONEDIMENSIONOPTIMIZATIONTOOLS_H_
