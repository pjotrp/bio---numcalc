//
// File: TwoPointsNumericalDerivative.cpp
// Created by: Julien Dutheil
// Created on: Mon May 28 10:33 2007
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "TwoPointsNumericalDerivative.h"

using namespace bpp;

void TwoPointsNumericalDerivative::updateDerivatives(const ParameterList & parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  if(_computeD1 && _variables.size() > 0)
  {
    if(_function1) _function1->enableFirstOrderDerivatives(false);
    if(_function2) _function2->enableSecondOrderDerivatives(false);
    _function->setParameters(parameters);
    _f1 = _function->getValue();
    string lastVar;
    for(unsigned int i = 0; i < _variables.size(); i++)
    {
      string var = _variables[i];
      if(parameters.getParameter(var) == NULL) continue;
      ParameterList p;
      if(i > 0)
      {
        vector<string> vars(2);
        vars[0] = var;
        vars[1] = lastVar;
        lastVar = var;
        p = parameters.subList(vars);
      }
      else
      {
        p = parameters.subList(var);
        lastVar = var;
      }
      double value = _function->getParameterValue(var);
      double h = (1 + std::abs(value)) * _h; 
      //Compute one other point:
      try
      {
        p[0]->setValue(value + h);
        _function->setParameters(p);
        _f2 = _function->getValue();
      }
      catch(ConstraintException & ce)
      {
        //Right limit raised, use backward approximation:
        try
        {
          p[0]->setValue(value - h);
          _function->setParameters(p);
          _f2 = _function->getValue();
          _der1[i] = (_f1 - _f2) / h;
        }
        catch(ConstraintException & ce)
        {
          //PB: can't compute derivative, because of a two narrow interval (lower than h)
          throw ce;
        }
      }
      //No limit raised, use forward approximation:
      _der1[i] = (_f2 - _f1) / h;
    }
    //Reset last parameter and compute analytical derivatives if any:
    if(_function1) _function1->enableFirstOrderDerivatives(_computeD1);
    _function->setParameters(parameters.subList(lastVar));
  }
  else
  {
    //Reset initial value and compute analytical derivatives if any.
    if(_function1) _function1->enableFirstOrderDerivatives(_computeD1);
    if(_function2) _function2->enableSecondOrderDerivatives(_computeD2);
    //Just in case derivatives are not computed:
    _function->setParameters(parameters);
    _f1 = _function->getValue();
  }
}

