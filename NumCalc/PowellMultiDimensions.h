//
// File: PowellMultiDimensions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov 17 15:16:45 2003
//

#ifndef _POWELLMULTIDIMENSIONS_H_
#define _POWELLMULTIDIMENSIONS_H_

#include "AbstractOptimizer.h"
#include "VectorTools.h"

class PowellMultiDimensions: public AbstractOptimizer
{
	/**************************************************************************/
	
	public:
		class PMDStopCondition: public AbstractOptimizationStopCondition
		{
			public:
				PMDStopCondition(PowellMultiDimensions * pmd);
				~PMDStopCondition();
			
			public:
				bool isToleranceReached() const;
		};
	
	/**************************************************************************/
	
	friend class PMDStopCondition;
		
	/**************************************************************************/

	protected:		
		class DirectionFunction: public Function
		{
			protected:
				mutable ParameterList _params, p;
				Vdouble xi;
				const Function * function;
			
			public:
				DirectionFunction(const Function * function);
				virtual ~DirectionFunction();
			
			public: // Function interface implementation:
				double f(const ParameterList & params) const;
				ParameterList getParameters() const;
			
			public: // Specific methods:
				void set(const ParameterList & p, const Vdouble & xi); 
		};
	
	/**************************************************************************/
	
	protected:
		double _fp;
		double _fret;
		ParameterList _pt;
		VVdouble _xi;
		
		int ncom;
		ParameterList pcom, xicom;
		DirectionFunction * f1dim;
		
	public:
		PowellMultiDimensions(const Function * function);
		virtual ~PowellMultiDimensions();
	
	public: // The Optimizer interface:
		
		void init(const ParameterList & params) throw (Exception); //redefinition
	
		double step() throw (Exception);
	
		double optimize() throw (Exception);

		double getFunctionValue() const;
	
	protected: // Specific function:
		
		void linmin(Vdouble & xi);
};


#endif	//_POWELLMULTIDIMENSIONS_H_
