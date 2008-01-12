//
// File: PowellMultiDimensions.h
// Created by: Julien Dutheil
// Created on: Mon Nov 17 15:16:45 2003
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

#ifndef _POWELLMULTIDIMENSIONS_H_
#define _POWELLMULTIDIMENSIONS_H_

#include "AbstractOptimizer.h"
#include "VectorTools.h"
#include "DirectionFunction.h"

namespace bpp
{

/**
 * @brief Powell's multi-dimensions optimization algorithm for one parameter.
 *
 * The code is an adaptation of the one discribed in:
 * <pre>
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * (ISBN 0-521-43108-5)
 * </pre>
 */
class PowellMultiDimensions:
  public AbstractOptimizer
{
	public:
		class PMDStopCondition:
      public AbstractOptimizationStopCondition
		{
			public:
				PMDStopCondition(PowellMultiDimensions * pmd):
          AbstractOptimizationStopCondition(pmd) {}
				virtual ~PMDStopCondition() {}

        PMDStopCondition * clone() const { return new PMDStopCondition(*this); }
			
			public:
				void init() {}
				bool isToleranceReached() const;
		};
	
	friend class PMDStopCondition;
		
	protected:
		double _fp;
		double _fret;
		ParameterList _pt;
		VVdouble _xi;
		
		unsigned int _ncom;
		ParameterList _pcom, _xicom;
		DirectionFunction _f1dim;
		
	public:
		PowellMultiDimensions(Function * function);
		virtual ~PowellMultiDimensions() {}

    PowellMultiDimensions * clone() const { return new PowellMultiDimensions(*this); }
	
	public:		
		
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */		
		double optimize() throw (Exception);
		/** @} */
		void doInit(const ParameterList & params) throw (Exception);
		
    double doStep() throw (Exception);	
	
};

} //end of namespace bpp.

#endif	//_POWELLMULTIDIMENSIONS_H_

