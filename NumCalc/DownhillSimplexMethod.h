//
// File: DownhillSimplexMethod.h
// Created by: Julien Dutheil
// Created on: Tue Nov  4 17:10:05 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _DOWNHILLSIMPLEXMETHOD_H_
#define _DOWNHILLSIMPLEXMETHOD_H_

#include "AbstractOptimizer.h"
#include "VectorTools.h"

// From the STL:
#include <cmath>
using namespace std;

namespace bpp
{

/**
 * @brief This implements the Downhill Simplex method in multidimensions.
 *
 * The code is an adaptation of the one discribed in:
 * <pre>
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * (ISBN 0-521-43108-5)
 * </pre>
 */
class DownhillSimplexMethod:
  public AbstractOptimizer
{
  public:
    class DSMStopCondition: public AbstractOptimizationStopCondition
    {
      public:
        DSMStopCondition(DownhillSimplexMethod * dsm):
          AbstractOptimizationStopCondition(dsm) {}
        virtual ~DSMStopCondition() {}

#ifndef NO_VIRTUAL_COV
        DSMStopCondition*
#else
        Clonable*
#endif
        clone() const { return new DSMStopCondition(*this); }
      
      public:
        void init() {}
        bool isToleranceReached() const;
    };
  
  friend class DSMStopCondition;
  
  protected:
    class Simplex:
      public vector<ParameterList>
  {
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
    unsigned int _iHighest, _iNextHighest, _iLowest;
  
  public:

    /**
     * @brief Build a new Downhill Simplex optimizer.
     *
     * @param function A pointer toward an object implementing the Optimizable interface.
     */
    DownhillSimplexMethod(Function * function);
  
    virtual ~DownhillSimplexMethod() {}

#ifndef NO_VIRTUAL_COV
    DownhillSimplexMethod*
#else 
    Clonable*
#endif
    clone() const { return new DownhillSimplexMethod(*this); }
  
  public:    
    /**
     * @name The Optimizer interface.
     *
     * @{
     */
    
    /**
     * @brief Multidimensional minimization of the function _function by the
     * downhill simplex method of Nelder and Mead.
     *
     * The simplex \f$p_{1..nDim+1,1..nDim}\f$.
     * is input. Its \f$nDim+1\f$ rows are nDim-dimensional vectors which are the vertices
     * of the starting simplex.
     * Also input is the vector \f$y_{1..nDim+1}\f$, whose components
     * must be preinitialized to the values of _function evaluated at the \f$nDim + 1\f$
     * vertices (rows).
     * On output, \f$p\f$ and \f$y\f$ will have been reset to \f$nDim + 1\f$
     * new points all within @c fTol of a minimum function value.
     */
    double optimize() throw (Exception);
    /** @} */

    void doInit(const ParameterList & params) throw (Exception);
    
    double doStep() throw (Exception);
  
  protected:
    
    /**
     * @name Specific inner methods
     *
     * @{
     */
    
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

    /** @} */
};

} //end of namespace bpp.

#endif  //_DOWNHILLSIMPLEXMETHOD_H_

