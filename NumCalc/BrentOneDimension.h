//
// File: BrentOneDimension.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 17 11:45:58 2003
//

#ifndef _BRENTONEDIMENSION_H_
#define _BRENTONEDIMENSION_H_

#include "AbstractOptimizer.h"

class BrentOneDimension: public AbstractOptimizer
{
	/**************************************************************************/
	
	public:
		class BODStopCondition: public AbstractOptimizationStopCondition
		{
			public:
				BODStopCondition(BrentOneDimension * bod);
				~BODStopCondition();
			
			public:
				bool isToleranceReached() const;
		};
	
	/**************************************************************************/
	
	friend class BODStopCondition;
	
	protected: // Fields:
		double a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
		double _xinf, _xsup;
		double _tolerance;

	public:
		BrentOneDimension(const Function * function);
		virtual ~BrentOneDimension();
	
	public: // Optimizer interface implemented here:
		
		/**
		 * @brief Initialize optimizer.
		 *
		 * The golden section search needs 2 initial guesses, so you must call the
		 * setInitialInterval() method first. This function actually performs:
		 * - Parameter list actualisation;
		 * - Initial bracketting;
		 * - Function evaluation count reseting.
		 *
		 * @param params The initial parameter list.
		 */
		void init(const ParameterList & params) throw (Exception); //redefinition
		double step() throw (Exception);
		double optimize() throw (Exception);
		double getFunctionValue() const;
	
	public: // Specific methods:
		
		void setInitialInterval(double inf, double sup);
		void setTolerance(double tolerance);
	
	public:
		
	static double CGOLD;
	static double ZEPS;

};


#endif	//_BRENTONEDIMENSION_H_
