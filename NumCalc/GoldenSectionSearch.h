//
// File: GoldenSectionSearch.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 10 10:42:17 2003
//

#ifndef _GOLDENSECTIONSEARCH_H_
#define _GOLDENSECTIONSEARCH_H_

#include "AbstractOptimizer.h"

class GoldenSectionSearch : public AbstractOptimizer
{	
	/**************************************************************************/
	
	public:
		class GSSStopCondition: public AbstractOptimizationStopCondition
		{
			public:
				GSSStopCondition(GoldenSectionSearch * gss);
				~GSSStopCondition();
			
			public:
				bool isToleranceReached() const;
		};
	
	/**************************************************************************/
	
	friend class GSSStopCondition;

	protected: // Fields:
		double f1, f2, x0, x1, x2, x3;
		double _xinf, _xsup;
	
	public: // Class constructor and destructor:
		
		GoldenSectionSearch(const Function * function);
		virtual ~GoldenSectionSearch();
	
	public: // Optimizer interface implemented here:
		
		/**
		 * @brief Initialize optimizer.
		 *
		 * The golden section search needs 2 initial guesses, so you must call the
		 * setInitialInterval() method first. this function actually performs:
		 * <ul>
		 * <li>Parameter list actualisation;</li>
		 * <li>Initial bracketting;</li>
		 * <li>Function evaluation count reseting.</li>
		 * </ul>
		 *
		 * @param params The initial parameter list.
		 */
		void init(const ParameterList & params) throw (Exception); //redefinition
		double step() throw (Exception);
		double optimize() throw (Exception);
		double getFunctionValue() const;
	
	public: // Specific methods:
		
		/**
		 * @brief Set the maximum number of iteration to perform during optimization.
		 *
		 * @param max The maximum number of itarations to perform.
		 */
		void setMaximumNumberOfIterations(int max);
	
		/**
		 * @brief Set the tolerance parameter. See the NRC book for more detail.
		 *
		 * The algorithm will stop when the size of the simplex is less than the tolerance parameter.
		 *
		 * @param tol The tolerance parameter.
		 */	
		void setTolerance(double tol);
	
		void setInitialInterval(double inf, double sup);
				
	protected: // Need these methods:
		
		static double R;
		static double C;
};


#endif	//_GOLDENSECTIONSEARCH_H_
