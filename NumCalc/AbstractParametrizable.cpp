//
// File: AbstractParametrizable.cpp
// Created by: Julien Dutheil
// Created on: Sun Mar 29 09:10 2009
//

/*
Copyright or © or Copr. CNRS, (November 19, 2004)

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

#include "AbstractParametrizable.h"

using namespace bpp;

AbstractParametrizable::AbstractParametrizable(const AbstractParametrizable & ap):
  _independentParameters(ap._independentParameters),
  _parameters(ap._parameters),
  _aliasListenersRegister(ap._aliasListenersRegister)
{
  //Actualize the register with adequate pointers:
  for(map<string, AliasParameterListener *>::iterator it = _aliasListenersRegister.begin();
    it != _aliasListenersRegister.end();
    it++)
  {
    it->second->setParameterList(&_parameters);
  }
}
   
AbstractParametrizable & AbstractParametrizable::operator=(const AbstractParametrizable & ap)
{
  _independentParameters = ap._independentParameters;
  _parameters = ap._parameters;
  _aliasListenersRegister = ap._aliasListenersRegister;

  //Actualize the register with adequate pointers:
  for(map<string, AliasParameterListener *>::iterator it = _aliasListenersRegister.begin();
    it != _aliasListenersRegister.end();
    it++)
  {
    it->second->setParameterList(&_parameters);
  }
  return *this;
}

AbstractParametrizable::~AbstractParametrizable()
{
  //Delete the registry content:
  for(map<string, AliasParameterListener *>::iterator it = _aliasListenersRegister.begin();
      it != _aliasListenersRegister.end();
      it++)
    delete it->second;
}

void AbstractParametrizable::aliasParameters(const string & p1, const string & p2)
  throw (ParameterNotFoundException, Exception)
{
  //In case this is the first time we call this method:
  if(_parameters.size() > 0 && _independentParameters.size() == 0) _independentParameters = _parameters;
 
  if(!_parameters.hasParameter(p1))
    throw ParameterNotFoundException("AbstractParametrizable::aliasParameters", p1);
  if(!_parameters.hasParameter(p2))
    throw ParameterNotFoundException("AbstractParametrizable::aliasParameters", p2);
  if(!_independentParameters.hasParameter(p2))
    throw Exception("AbstractParametrizable::aliasParameters. Parameter " + p2 + " is already aliased to a parameter and can't be aliased twice.");

  string id = "__alias_" + p2 + "_to_" + p1;
  string idCheck = "__alias_" + p1 + "_to_" + p2;
  if(_aliasListenersRegister.find(idCheck) != _aliasListenersRegister.end())
    throw Exception("AbstractParametrizable::aliasParameters. Trying to alias parameter " + p2 + " to " + p1 + ", but parameter " + p1 + " is already aliased to parameter " + p2 + ".");
  Parameter *param1 = _parameters.getParameter(p1);
  Parameter *param2 = _parameters.getParameter(p2);
  if(!param1->hasConstraint())
  {
    if(param2->hasConstraint())
      throw Exception("AbstractParametrizable::aliasParameters. Cannot alias parameter " + p2 + " to " + p1 + ", because the constraints aatached to these two parameters are different.");
  }
  else
    //We use a small trick here, we test the constraints on the basis of their string description (C++ does not provide a default operator==() :( ).
    if(param2->hasConstraint() && (param1->getConstraint()->getDescription() != param2->getConstraint()->getDescription()))
        throw Exception("AbstractParametrizable::aliasParameters. Cannot alias parameter " + p2 + " to " + p1 + ", because the constraints aatached to these two parameters are different.");

  //Every thing seems ok, let's create the listener and register it:
  AliasParameterListener * aliasListener = new AliasParameterListener(id, _parameters.whichParameterHasName(p2), &_parameters);
  _aliasListenersRegister[id] = aliasListener;
  //Now we add it to the appropriate parameter, that is p1.
  //The parameter will not own the listener, the bookkeeping being achieved by the register:
  param1->addParameterListener(aliasListener, false);
  //Finally we remove p2 from the list of independent parameters:
  _independentParameters.deleteParameter(p2);
}

void AbstractParametrizable::unaliasParameters(const string & p1, const string & p2)
  throw (ParameterNotFoundException, Exception)
{
  if(!_parameters.hasParameter(p1))
    throw ParameterNotFoundException("AbstractParametrizable::unaliasParameters", p1);
  if(!_parameters.hasParameter(p2))
    throw ParameterNotFoundException("AbstractParametrizable::unaliasParameters", p2);

  string id = "__alias_" + p2 + "_to_" + p1;
  map<string, AliasParameterListener *>::iterator it = _aliasListenersRegister.find(id);
  if(it == _aliasListenersRegister.end())
    throw Exception("AbstractParametrizable::unaliasParameters. Parameter " + p2 + " is not aliased to parameter " + p1 + ".");
  //Remove the listener:
  _parameters.getParameter(p1)->removeParameterListener(id);
  delete it->second;
  _aliasListenersRegister.erase(it);
  //Finally we re-add p2 to the list of independent parameters:
  _independentParameters.addParameter(*_parameters.getParameter(p2));
}
	
