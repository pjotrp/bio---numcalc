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
				mutable ParameterList _params, p, _xt;
				Vdouble xi;
				Function * function;
			
			public:
				DirectionFunction(Function * function);
				virtual ~DirectionFunction();
			
			public: // Function interface implementation:
				void setParameters(const ParameterList & params)
					throw (ParameterNotFoundException, ConstraintException);
				double getValue() const throw (Exception);
				ParameterList getParameters() const throw (Exception);
				double getParameter(const string & name) const throw (ParameterNotFoundException) { return 0; };
				void setAllParametersValues(const ParameterList & params) 
					throw (ParameterNotFoundException, ConstraintException) {}
				void setParameterValue(const string & name, double value) 
					throw (ParameterNotFoundException, ConstraintException) {}
				void setParametersValues(const ParameterList & params)
					throw (ParameterNotFoundException, ConstraintException) {}
				void matchParametersValues(const ParameterList & params)
					throw (ConstraintException) {};
			
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
		PowellMultiDimensions(Function * function);
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
