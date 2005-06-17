//
// File: Optimizer.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov  4 16:01:27 2003
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

#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "Functions.h"
#include "ParameterList.h"
#include "OptimizationStopCondition.h"

// From the STL:
#include <iostream>

using namespace std;

/**
 * @brief This is the base interface for all optimization methods.
 * 
 * An optimizer will deal with object that implement the Optimizable interface.
 */
class Optimizer
{
	public:
		Optimizer() {}
		virtual ~Optimizer() {}
	
	public:
		
		/**
		 * @brief Set the initial values of the parameters.
		 *
		 * @param params The initial values of parameters.
		 */
		virtual void init(const ParameterList & params) throw (Exception) = 0;

		/**
		 * @brief Perform a step optimization.
		 *
		 * @return the value of the function after this step.
		 */
		virtual double step() throw (Exception) = 0;

		/**
		 * @brief Get the current parameters.
		 *
		 * @return The current parameters.
		 */
		virtual ParameterList getParameters() const = 0;
		
		/**
		 * @brief Get the corresponding function evaluation.
		 *
		 * @return The value of the function at the point specified by _parameters.
		 */
		virtual double getFunctionValue() const = 0;
		
		/**
		 * @brief Perform as many optimization steps untill the stop condition is met.
		 *
		 * @return The value of the function after optimization is completed.
		 */
		virtual double optimize() throw (Exception) = 0;
	
		/**
		 * @brief Set the function to optimize.
		 *
		 * @param function The function to optimize.
		 */
		virtual void setFunction(Function * function) = 0;
		
		/**
		 * @brief Get the whole Optimizable object.
		 *
		 * @return A const pointer toward the object to be optimized.
		 */
		virtual const Function * getFunction() const = 0;

		/**
		 * @brief Set the message handler for this optimizer.
		 *
		 * The message handler keeps all messages that the optimizer may send.
		 * The default handler is set to standard output, but you can pass any
		 * ostream object (cerr, ofstream, etc.).
		 *
		 * @param mh The message handler to use.
		 */
		virtual void setMessageHandler(ostream * mh) = 0;
		
		/**
		 * @brief Set the profiler for this optimizer.
		 *
		 * The profiler keeps all the intermediate values taken by the parameters.
		 * The default profiler is set to standard output, but you can pass any
		 * ostream object (cerr, ofstream, etc.).
		 *
		 * @param profiler The profiler to use.
		 */
		virtual void setProfiler(ostream * profiler) = 0;
		
		/**
		 * @brief Get the number of function evaluation done since call of
		 * the init function.
		 *
		 * @return the number of function evaluations.
		 */
		virtual	int getNumberOfEvaluations() const = 0;
		
		/**
		 * @brief Set the stopping condition of the optimization algorithm.
		 *
		 * @param stopCondition The stoping condition to use while optimizing.
		 */
		virtual void setStopCondition(OptimizationStopCondition * stopCondition) = 0;

		/**
		 * @brief Get the stopping condition of the optimization algorithm.
		 *
		 * @return The stopping condition used while optimizing.
		 */
		virtual OptimizationStopCondition * getStopCondition() = 0;

		/**
		 * @brief Get the default stopping condition of the optimization algorithm.
		 *
		 * @return The default stopping condition used while optimizing.
		 */
		virtual OptimizationStopCondition * getDefaultStopCondition() = 0;
		
		/**
		 * @brief Tell if the tolerance level is reached.
		 *
		 * @return Whether the tolerance is reached or not.
		 * @see OptimizationStopCondition
		 */
		virtual bool isToleranceReached() const = 0;
		
		/**
		 * @brief Tell if the maximum number of function evaluations is reached.
		 *
		 * @return Whether the maximum number of function evaluations is reached or not.
		 */
		virtual bool isMaximumNumberOfEvaluationsReached() const = 0;

		/**
		 * @brief Set the verbose level.
		 *
		 * 0 = off
		 * 1 = on
		 * 2 = more verbose
		 * 3 = even more, etc.
		 *
		 * In most cases, only the 0 and 1 levels are implemented.
		 *
		 * @param v verbose level.
		 */
		virtual void setVerbose(unsigned int v) = 0;

		/**
		 * @brief Get the verbose level.
		 *
		 * @return verbose level.
		 */
		virtual unsigned int getVerbose() const = 0;
};

#endif	//_OPTIMIZER_H_
