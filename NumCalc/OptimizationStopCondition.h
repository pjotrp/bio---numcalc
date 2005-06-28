//
// File: OptimizationStopCondition.h
// Created by: Julien Dutheil
// Created on: Tue Dec 23 11:51:31 2003
//

/*
Copyright ou © ou Copr. CNRS, (19 Novembre 2004) 

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
Copyright or © or Copr. CNRS, (November 19, 2004)

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

#ifndef _OPTIMIZATIONSTOPCONDITION_H_
#define _OPTIMIZATIONSTOPCONDITION_H_

#include "ParameterList.h"

class Optimizer;
	
/******************************************************************************/

class OptimizationStopCondition
{
	public:
		OptimizationStopCondition();
		virtual ~OptimizationStopCondition();
	
	public:

		/**
		 * @brief Initialize the condition.
		 */
		virtual void init() = 0;

		/**
		 * @brief Tell if the we reached the desired tolerance with a given 
		 * new set of estimates.
		 *
		 * The new parameter list is compared to the last estimates,
		 * and the lastParameterEstimates list is actulaized with the newParameters list.
		 *
		 * @return True if the tolerance level is reached.
		 */
		virtual bool isToleranceReached() const = 0;
		
		/**
		 * @brief Set the tolerance parameter.
		 *
		 * @param tolerance The tolerance parameter.
		 */	
		virtual void setTolerance(double tolerance) = 0;

		/**
		 * @brief Get the tolerance parameter.
		 *
		 * @return The tolerance parameter.
		 */	
		virtual double getTolerance() const = 0;
};

/******************************************************************************/

class AbstractOptimizationStopCondition: public virtual OptimizationStopCondition
{
	protected:
		const Optimizer * _optimizer;
		double _tolerance;

		/**
		 * @brief Count the number of times the isToleranceReached() function
		 * has been called.
		 */
		mutable double _callCount;
	
		int _burnin;

	
	public:
		AbstractOptimizationStopCondition(const Optimizer * optimizer);
		AbstractOptimizationStopCondition(const Optimizer * optimizer, double tolerance);
		AbstractOptimizationStopCondition(const Optimizer * optimizer, int burnin);
		AbstractOptimizationStopCondition(const Optimizer * optimizer, double tolerance, int burnin);
	
		virtual ~AbstractOptimizationStopCondition();

	public:
		void setTolerance(double tolerance);
		double getTolerance() const;
	
		virtual void resetCounter();
		virtual void setBurnin(int burnin);
		virtual int getBurnin() const;

};
	
/******************************************************************************/

class ParametersStopCondition: public virtual AbstractOptimizationStopCondition
{
	protected:
		/**
		 * @brief The last estimates of the parameters.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable ParameterList _lastParametersEstimates;
		
		/**
		 * @brief The new estimates of the parameters.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable ParameterList _newParametersEstimates;
	
	public:
		ParametersStopCondition(const Optimizer * optimizer);
		ParametersStopCondition(const Optimizer * optimizer, double tolerance);
		ParametersStopCondition(const Optimizer * optimizer, int burnin);
		ParametersStopCondition(const Optimizer * optimizer, double tolerance, int burnin);
		
		~ParametersStopCondition();
	
	public:
		void init();

		bool isToleranceReached() const;
};

/******************************************************************************/

class FunctionStopCondition: public virtual AbstractOptimizationStopCondition
{
	protected:
		/**
		 * @brief The last value of the function.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable double _lastFunctionValue;
		
		/**
		 * @brief The new value of the function.
		 *
		 * This is used by the isToleranceReached() method.
		 */
		mutable double _newFunctionValue;
	
	public:
		FunctionStopCondition(const Optimizer * optimizer);
		FunctionStopCondition(const Optimizer * optimizer, double tolerance);
		FunctionStopCondition(const Optimizer * optimizer, int burnin);
		FunctionStopCondition(const Optimizer * optimizer, double tolerance, int burnin);
		
		~FunctionStopCondition();
	
	public:
		void init();
		bool isToleranceReached() const;

};

/******************************************************************************/

#endif	//_OPTIMIZATIONSTOPCONDITION_H_
