//
// File: Functions.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Nov  9 23:11:00 2003
//

/*
Copyright ou © ou Copr. CNRS, (17 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour le calcul numérique.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

Julien.Dutheil@univ-montp2.fr

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

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "ParameterList.h"
#include "Parametrizable.h"
#include "ParameterExceptions.h"

/**
 * @brief This is the function abstract class.
 */
class Function : public virtual Parametrizable
{		
	public:
		Function() {}
		virtual ~Function() {}

	public:

		/**
		 * @brief Set the point where the function must be computed.
		 *
		 * @param parameters The parameter set to pass to the function.
		 */
		virtual void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) = 0;

		/**
		 * @brief Get the current point at which the function is evaluated.
		 *
		 * @return The set of parameters corresponding to the current point.
		 * @throw Exception If no point is defined.
		 */
		virtual ParameterList getParameters() const throw (Exception) = 0;

		/**
		 * @brief Get the value of the function at the current point.
		 *
		 * @return The value of the function.
		 * @throw Exception If no point is specified or if an error occured.
		 */
		virtual double getValue() const throw (Exception) = 0;
		
		/**
		 * @brief Get the value of the function according to a given set of parameters.
		 * 
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 * @throw Exception If an error occured.
		 */
		virtual double f(const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getValue();
		}
};

/**
 * @brief This is the interface for first order derivable functions.
 */
class DerivableFirstOrder : public virtual Function
{
	public:
		DerivableFirstOrder() {}
		virtual ~DerivableFirstOrder() {}

	public:
		
		/**
		 * @brief Get the derivative of the function at the current point.
		 *
		 * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{df}{dx} @f$.
		 * @return The value of the function.
		 * @throw Exception If no point is specified or if an error occured.
		 */
		virtual double getFirstOrderDerivative(const string & variable) const throw (Exception) = 0;
		
		/**
		 * @brief Get the value of the first derivative of the function
		 * according to a given set of parameters.
		 *
		 * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{df}{dx} @f$.
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 * @throw Exception If an error occured.
		 */
		virtual double df(const string & variable, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getFirstOrderDerivative(variable);
		}
};

/**
 * @brief This is the interface for second order derivable functions.
 */
class DerivableSecondOrder : public virtual DerivableFirstOrder
{
	public:
		DerivableSecondOrder() {}
		virtual ~DerivableSecondOrder() {}

	public:

		/**
		 * @brief Get the second order derivative of the function at the current point.
		 *
		 * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x^2} @f$.
		 * @return The value of the function.
		 * @throw Exception If no point is specified or if an error occured.
		 */
		virtual double getSecondOrderDerivative(const string & variable) const throw (Exception) = 0;
	
		/**
		 * @brief Get the value of the second order derivative of the function
		 * according to a given set of parameters.
		 *
		 * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x^2} @f$.
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 * @throw Exception If an error occured.
		 */
		virtual double d2f(const string & variable, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getSecondOrderDerivative(variable);
		}		

		/**
		 * @brief Get the value of the cross derivative of the function
		 * according to a given set of parameters.
		 *
		 * @param variable1  The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
		 * @param variable2  The name of the @f$ y @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
		 * @return The value of the function with the given parameter set.
		 * @throw Exception If an error occured.
		 */
		virtual double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) = 0;	
		
		/**
		 * @brief Get the value of the cross derivative of the function
		 * according to a given set of parameters.
		 *
		 * @param variable1  The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
		 * @param variable2  The name of the @f$ y @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
		 * @param parameters The parameter set to pass to the function.
		 * @return The value of the function with the given parameter set.
		 * @throw Exception If an error occured.
		 */
		virtual double d2f(const string & variable1, const string & variable2, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getSecondOrderDerivative(variable1, variable2);
		}
};




class InfinityFunctionWrapper : public virtual Function
{
	protected:
		Function * _function;
		mutable bool _constraintMatch;
		
	public:
		InfinityFunctionWrapper(Function * function): _function(function), _constraintMatch(false) {}
		virtual ~InfinityFunctionWrapper() {}

	public:

		void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setParameters(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}

		ParameterList getParameters() const throw (Exception) {
			return _function -> getParameters();	
		}

		double getValue() const throw (Exception) {
			return _constraintMatch ? -log(0.) :	_function -> getValue();
		}
		
		double f(const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getValue();
		}
		
		double getParameter(const string & name) const throw (ParameterNotFoundException) {
			return _function -> getParameter(name);
		}
			
		void setAllParametersValues(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setAllParametersValues(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}
		
		void setParameterValue(const string & name, double value) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setParameterValue(name, value);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}
		
		void setParametersValues(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException) {
			try {
				_function -> setParametersValues(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}
		
		void matchParametersValues(const ParameterList & parameters) throw (ConstraintException) {
			try {
				_function -> matchParametersValues(parameters);
				_constraintMatch = false;
			} catch(ConstraintException ce) {
				_constraintMatch = true;
			}
		}

};


class InfinityDerivableFirstOrderWrapper : public virtual InfinityFunctionWrapper
{
	public:
		InfinityDerivableFirstOrderWrapper(DerivableFirstOrder * function): InfinityFunctionWrapper(function) {}
		virtual ~InfinityDerivableFirstOrderWrapper() {}

	public:
		
		double getFirstOrderDerivative(const string & variable) const throw (Exception) {
			return _constraintMatch ? -log(0.) :	(dynamic_cast<DerivableFirstOrder *>(_function) -> getFirstOrderDerivative(variable));		
		}
		
		double df(const string & variable, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getFirstOrderDerivative(variable);
		}
};


class InfinityDerivableSecondOrderWrapper : public virtual InfinityDerivableFirstOrderWrapper
{
	public:
		InfinityDerivableSecondOrderWrapper(DerivableFirstOrder * function):
			InfinityFunctionWrapper(function),
			InfinityDerivableFirstOrderWrapper(function) {}
		virtual ~InfinityDerivableSecondOrderWrapper() {}

	public:

		double getSecondOrderDerivative(const string & variable) const throw (Exception) {
			return _constraintMatch ? -log(0.) :	(dynamic_cast<DerivableSecondOrder *>(_function) -> getSecondOrderDerivative(variable));					
		}
	
		double d2f(const string & variable, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getSecondOrderDerivative(variable);
		}		

	  double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) {
			return _constraintMatch ? -log(0.) :	(dynamic_cast<DerivableSecondOrder *>(_function) -> getSecondOrderDerivative(variable1, variable2));			
		}
		
		double d2f(const string & variable1, const string & variable2, const ParameterList & parameters) throw (Exception) {
			setParameters(parameters);
			return getSecondOrderDerivative(variable1, variable2);
		}

};

#endif	//_FUNCTIONS_H_
