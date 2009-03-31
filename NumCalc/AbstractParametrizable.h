//
// File: AbstractParametrizable.h
// Created by: Julien Dutheil
// Created on: Sun Mar 29 09:10 2009
// Created from file Parametrizable.h
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

#ifndef _ABSTRACTPARAMETRIZABLE_H_
#define _ABSTRACTPARAMETRIZABLE_H_

#include "Parametrizable.h"

//From the STL:
#include <map>
using namespace std;

namespace bpp
{

/**
 * @brief A partial implementation of the Parametrizable interface.
 *
 * Parameters are stored in a protected ParameterList object.
 *
 * The abstract fireParameterChanged() method is provided so that the derived class
 * know when a parameter has changed, and can be updated.
 * All methods call the corresponding method in ParameterList and then call the
 * fireParameterChanged() method.
 */
class AbstractParametrizable:
  public virtual Parametrizable
{
  private:
    /**
     * This lst is only for bookkeeping only.
     * It will automatically be updated from _parameters.
     */
		mutable ParameterList _independentParameters;

	protected:
		mutable ParameterList _parameters;

    class AliasParameterListener:
      public ParameterListener
    {
      protected:
        string _id;
        unsigned int _alias;
        ParameterList *_pl;
        string _name;

      public:
        AliasParameterListener(const string & id, unsigned int alias, ParameterList * pl): _id(id), _alias(alias), _pl(pl)
        {
          //This allow us to check if the parameter position have changed at some point...
          _name = (*_pl)[alias]->getName();
        }
        AliasParameterListener * clone() const { return new AliasParameterListener(*this); }

      public:
        const string & getId() const { return _id; }

        void setParameterList(ParameterList *pl) { _pl = pl; }

        void parameterNameChanged(ParameterEvent & event) throw (Exception)
        {
          Parameter * p = (*_pl)[_alias];
          if(p->getName() != _name)
            throw Exception("AbstractParametrizable::AliasParameterListener::parameterNameChanged. Error, aliased parameter have change, maybe because a parameter was removed?");
          p->setName(event.getParameter()->getName());
          _name = p->getName();
        }
    
        void parameterValueChanged(ParameterEvent & event) throw (Exception)
        {
          Parameter * p = (*_pl)[_alias];
          if(p->getName() != _name)
            throw Exception("AbstractParametrizable::AliasParameterListener::parameterValueChanged. Error, aliased parameter have change, maybe because a parameter was removed?");
          p->setValue(event.getParameter()->getValue());
        }
      
    };

    /**
     * Contains all parameter listeners for maintening alias relationships.
     * The registry will be updated appropriately upon cloning and deleting.
     */
    map<string, AliasParameterListener *> _aliasListenersRegister;
	
	public:
		AbstractParametrizable() {}

    AbstractParametrizable(const AbstractParametrizable & ap);
    
    AbstractParametrizable & operator=(const AbstractParametrizable & ap);

		virtual ~AbstractParametrizable();

	public:

		const ParameterList & getParameters() const { return _parameters; }
		
    const ParameterList & getIndependentParameters() const
    {
      //A small precaution in case _independentParameters was not initialized:
      if(_parameters.size() > 0 && _independentParameters.size() == 0) _independentParameters = _parameters;
      return _independentParameters;
    }
    
    const Parameter & getParameter(const string & name) const throw (ParameterNotFoundException)
    {
      const Parameter * p = _parameters.getParameter(name);
      if(p) return *p;
      else throw ParameterNotFoundException("AbstractParametrizable::getParameter.", name);
    }
	
		double getParameterValue(const string & name) const
			throw (ParameterNotFoundException)
		{ 
			return _parameters.getParameter(name)->getValue();
		}

		void setAllParametersValues(const ParameterList & parameters) 
			throw (ParameterNotFoundException, ConstraintException)
		{
			_parameters.setAllParametersValues(parameters);
      _independentParameters.matchParametersValues(_parameters);
			fireParameterChanged(parameters);
		}

		void setParameterValue(const string & name, double value) 
			throw (ParameterNotFoundException, ConstraintException)
		{ 
			_parameters.setParameterValue(name, value);
      _independentParameters.matchParametersValues(_parameters);
			fireParameterChanged(_parameters.subList(name));
		}

		void setParametersValues(const ParameterList & parameters)
			throw (ParameterNotFoundException, ConstraintException)
		{ 
			_parameters.setParametersValues(parameters);
      _independentParameters.matchParametersValues(_parameters);
			fireParameterChanged(parameters);
		}

		void matchParametersValues(const ParameterList & parameters)
			throw (ConstraintException)
		{ 
			_parameters.matchParametersValues(parameters);
      _independentParameters.matchParametersValues(_parameters);
			fireParameterChanged(parameters);
		}

    unsigned int getNumberOfParameters() const { return _parameters.size(); }
    
    unsigned int getNumberOfIndependentParameters() const { return _independentParameters.size(); }

    void aliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception);

    void unaliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception);
	
		/**
		 * @brief Notify the class when one or several parameters have changed.
		 *
		 * @param parameters A ParameterList object with parameters that changed.
		 */
		virtual void fireParameterChanged(const ParameterList & parameters) = 0;

};

} //end of namespace bpp.
#endif //_ABSTRACTPARAMETRIZABLE_H_

