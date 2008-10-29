//
// File: ThreePointsNumericalDerivative.cpp
// Created by: Julien Dutheil
// Created on: Thu Aug 17 15:00 2006
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

#include "ThreePointsNumericalDerivative.h"

using namespace bpp;

void ThreePointsNumericalDerivative::setParameters(const ParameterList & parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  if(_computeD1 && _variables.size() > 0)
  {
    if(_function1) _function1->enableFirstOrderDerivatives(false);
    if(_function2) _function2->enableSecondOrderDerivatives(false);
    _function->setParameters(parameters);
    _f2 = _function->getValue();
    string lastVar;
    ParameterList p;
    for(unsigned int i = 0; i < _variables.size(); i++)
    {
      string var = _variables[i];
      if(parameters.getParameter(var) == NULL) continue;
      if(i > 0)
      {
        vector<string> vars(2);
        vars[0] = var;
        vars[1] = lastVar;
        p = parameters.subList(vars);
      }
      else
      {
        p = parameters.subList(var);
      }
      lastVar = var;
      double value = _function->getParameterValue(var);
      double h = (1. + std::abs(value)) * _h; 
      //Compute one other point:
      try
      {
        p[0]->setValue(value - h);
        _function->setParameters(p); //also reset previous parameter...
        p = p.subList(0);
        _f1 = _function->getValue();
        try
        {
          p[0]->setValue(value + h);
          _function->setParameters(p);
          _f3 = _function->getValue();
        }
        catch(ConstraintException & ce)
        {
          //Right limit raised, use backward approximation:
          p[0]->setValue(value - h);
          _function->setParameters(p);
          _f1 = _function->getValue();
          p[0]->setValue(value - 2*h);
          _function->setParameters(p);
          _f3 = _function->getValue();
          _der1[i] = (_f2 - _f1) / h;
          _der2[i] = (_f2 - 2.*_f1 + _f3) / (h*h);        
        }
        //No limit raised, use central approximation:
        _der1[i] = (-_f1 + _f3) / (2.*h);
        _der2[i] = (_f1 -2*_f2 + _f3) / (h*h);
      }
      catch(ConstraintException & ce)
      {
        //Left limit raised, use forward approximation:
        p[0]->setValue(value + h);
        _function->setParameters(p);
        _f3 = _function->getValue();
        p[0]->setValue(value + 2*h);
        _function->setParameters(p);
        _f1 = _function->getValue();
        _der1[i] = (_f3 - _f2) / h;
        _der2[i] = (_f1 - 2.*_f3 + _f2) / (h*h);
      }
    }

    if(_computeCrossD2)
    {
      string lastVar1, lastVar2;
      for(unsigned int i = 0; i < _variables.size(); i++)
      {
        string var1 = _variables[i];
        if(parameters.getParameter(var1) == NULL) continue;
        for(unsigned int j = 0; j < _variables.size(); j++)
        {
          if(j==i)
          {
            _crossDer2(i,j) = _der2[i];
            continue;
          }
          string var2 = _variables[j];
          if(parameters.getParameter(var2) == NULL) continue;
        
          vector<string> vars(2);
          vars[0] = var1;
          vars[1] = var2;
          if(i > 0 && j > 0)
          {
            if(lastVar1 != var1 && lastVar1 != var2) vars.push_back(lastVar1);
            if(lastVar2 != var1 && lastVar2 != var2) vars.push_back(lastVar2);
          }
          p = parameters.subList(vars);
        
          double value1 = _function->getParameterValue(var1);
          double value2 = _function->getParameterValue(var2);
          double h1 = (1. + std::abs(value1)) * _h; 
          double h2 = (1. + std::abs(value2)) * _h; 
        
          //Compute 4 additional points:
          try
          {
            p[0]->setValue(value1 - h1);
            p[1]->setValue(value2 - h2);
            _function->setParameters(p); //also reset previous parameter...
            vector<unsigned int> tmp(2);
            tmp[0] = 0;
            tmp[1] = 1;
            p = p.subList(tmp); //removed the previous parameters.
            _f11 = _function->getValue();

            p[1]->setValue(value2 + h2);
            _function->setParameters(p.subList(1));
            _f12 = _function->getValue();

            p[0]->setValue(value1 + h1);
            _function->setParameters(p.subList(0));
            _f22 = _function->getValue();

            p[1]->setValue(value2 - h2);
            _function->setParameters(p.subList(1));
            _f21 = _function->getValue();

            _crossDer2(i,j) = ((_f22 - _f21) - (_f12 - _f11)) / (4*h1*h2);
          }
          catch(ConstraintException & ce)
          {
            throw Exception("ThreePointsNumericalDerivative::setParameters. Could not compute cross derivatives at limit.");
          }

          lastVar1 = var1;
          lastVar2 = var2;
        }
      }
    }
   
    //Reset last parameter and compute analytical derivatives if any.
    if(_function1) _function1->enableFirstOrderDerivatives(_computeD1);
    if(_function2) _function2->enableSecondOrderDerivatives(_computeD2);
    _function->setParameters(parameters.subList(lastVar));
  }
  else
  {
    //Reset initial value and compute analytical derivatives if any.
    if(_function1) _function1->enableFirstOrderDerivatives(_computeD1);
    if(_function2) _function2->enableSecondOrderDerivatives(_computeD2);
    _function->setParameters(parameters);
    //Just in case derivatives are not computed:
    _f2 = _function->getValue();
  }
}

