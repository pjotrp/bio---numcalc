//
// File: DownhillSimplexMethod.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov  4 17:10:05 2003
//

#ifndef _DOWNHILLSIMPLEXMETHOD_H_
#define _DOWNHILLSIMPLEXMETHOD_H_

#include "AbstractOptimizer.h"
#include "VectorTools.h"

// From the STL:
#include <cmath>

/**
 * @brief This implements the Downhill Simplex method in multidimensions.
 * The code is an adaptation of the one discribed in:
 * <pre>
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * (ISBN 0-521-43108-5)
 * </pre>
 */
class DownhillSimplexMethod: public AbstractOptimizer
{
	/**************************************************************************/
	
	public:
		class DSMStopCondition: public AbstractOptimizationStopCondition
		{
			public:
				DSMStopCondition(DownhillSimplexMethod * dsm);
				~DSMStopCondition();
			
			public:
				bool isToleranceReached() const;
		};
	
	/**************************************************************************/
	
	friend class DSMStopCondition;
	
	protected:
		class Simplex: public vector<ParameterList> {
			public: // Class constructor and destructor:
				Simplex();
				virtual ~Simplex();
			
			public: // Methods:
				virtual int getDimension() const;
		};
		
	protected:
		Simplex _simplex;
		Vdouble _y;
		ParameterList _pSum;
		int _iHighest, _iNextHighest, _iLowest;
	
	public: // constructor and destructor:

		/**
		 * @brief Build a new Downhill Simplex optimizer.
		 *
		 * @param function A pointer toward an object implementing the Optimizable interface.
		 */
		DownhillSimplexMethod(Function * function);
	
		virtual ~DownhillSimplexMethod();
	
	public: // The Optimizer interface:
		
		void init(const ParameterList & params) throw (Exception); //redefinition
	
		double step() throw (Exception);
	
		/**
		 * @brief Multidimensional minimization of the function _function by the
		 * downhill simplex method of Nelder and Mead.
		 *
		 * The simplex <i>p</i>[1..nDim+1][1..nDim]
		 * is input. Its <i>nDim+1</i> rows are nDim-dimensional vectors which are the vertices
		 * of the starting simplex. Also input is the vector <i>y</i>[1..nDim+1], whose components
		 * must be preinitialized to the values of _function evaluated at the <i>nDim + 1</i>
		 * vertices (rows). On output, <i>p</i> and <i>y</i> will have been reset to <i>nDim + 1</i>
		 * new points all within <i>fTol</i> of a minimum function value.
		 */
		double optimize() throw (Exception);
		double getFunctionValue() const;
	
	protected:
		
		/**
		 * @brief Update the _pSum variable.
		 */
		ParameterList getPSum();
	
		/**
		 * @brief Extrapolates by a factor fac throough the face of the simplex
		 * from the high point, try it, an dreplaces the high point if the new point is better.
		 *
		 * @param fac Extrapolation factor.
		 * @return The value of the function for the new point.
		 */
		double amotry(double fac);
};


#endif	//_DOWNHILLSIMPLEXMETHOD_H_
