//
// File: AbstractOptimizer.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Dec 22 12:18:09 2003
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

#ifndef _ABSTRACTOPTIMIZER_H_
#define _ABSTRACTOPTIMIZER_H_

#include "Optimizer.h"

/**
 * @brief This class offers a basic implementation of the Optimizer interface.
 */
class AbstractOptimizer: public virtual Optimizer
{
	protected:
		
		/**
		 * @brief The function to optimize.
		 */
		Function * _function;
	
		/**
		 * @brief The parameters that will be optimized.
		 */
		ParameterList _parameters;
	
		/**
		 * @brief The message handler.
		 */
		ostream * _messageHandler;
	
		/**
		 * @brief The profiler.
		 */
		ostream * _profiler;
		
		/**
		 * @brief The constraint policy.
		 *
		 * Must be one the following:
		 * - CONSTRAINTS_KEEP: keep the constraint associated to the parameters (default).
		 * - CONSTRAINTS_IGNORE: remove all constraints.
		 * - CONSTRAINTS_AUTO: use AutoParameters to deal with constraints.
		 *
		 * @see AutoParameter
		 */ 		 
		string _constraintPolicy;

		/**
		 * @brief Tell if the tolerance level has been reached.
		 *
		 * This field is initilaised by the init() method, maintained by the
		 * step() method and used in the optimize() method.
		 */
		bool _tolIsReached;
		
		/**
		 * @brief The stoping condition to use while optimizing.
		 */
		OptimizationStopCondition * _stopCondition;
		
		/**
		 * @brief The default stoping condition to use while optimizing.
		 */
		OptimizationStopCondition * _defaultStopCondition;

		/**
		 * @brief The maximum number of function evaluations allowed.
		 */
		int _nbEvalMax;
		
		/**
		 * @brief The current number of function evaluations achieved.
		 */
		int _nbEval;

		/**
		 * @brief State of the verbose mode: > 0 = enabled.
		 *
		 * This may not be used by the Optimizer.
		 */
		unsigned int _verbose;


	public:
		AbstractOptimizer(Function * function);
		virtual ~AbstractOptimizer();
	
	public:
		
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */
		void init(const ParameterList & params) throw (Exception);
		ParameterList getParameters() const;
		void setFunction(Function * function) { _function = function; if(function != NULL) _stopCondition -> init(); }
		const Function * getFunction() const { return _function; }
		double getFunctionValue() const throw (Exception)
		{
			if(_function == NULL) throw Exception("AbstractOptimizer::getFunctionValue. No function associated to this optimizer.");
			return _function -> getValue();
		}
		void setMessageHandler(ostream * mh);
		void setProfiler(ostream * profiler);
		int getNumberOfEvaluations() const;
		void setStopCondition(OptimizationStopCondition * stopCondition);
		OptimizationStopCondition * getStopCondition();
		OptimizationStopCondition * getDefaultStopCondition();
		bool isToleranceReached() const;
		bool isMaximumNumberOfEvaluationsReached() const;
		void setVerbose(unsigned int v) { _verbose = v; }
		unsigned int getVerbose() const { return _verbose; }
		/** @} */
	
	protected: // Some util:
		
		/**
		 * @brief Build a list of AutoParameter instead of Parameter.
		 */
		void autoParameter();
	
		/**
		 * @brief Remove the constraints of all the arguments.
		 */
		void ignoreConstraints();
	
		/**
		 * @brief Print to the profile if there is one.
		 *
		 * @param v The double value to print.
		 */
		void profile(double v);
	
		/**
		 * @brief Print to the profile if there is one.
		 *
		 * @param s The string to print to the profile.
		 */
		void profile(const string & s);
	
		/**
		 * @brief Print to the profile if there is one and end line.
		 *
		 * @param v The double value to print.
		 */
		void profileln(double v);
	
		/**
		 * @brief Print to the profile if there is one and end line.
		 *
		 * @param s The string to print to the profile.
		 */
		void profileln(const string & s);
	
		/**
		 * @brief Print parameters and corresponding function evaluation to profiler.
		 *
		 * @param params The parameters to print.
		 * @param value  The function evaluation.
		 */
		void printPoint(const ParameterList & params, double value);
		
		/**
		 * @brief Give a message to print to the message handler.
		 *
		 * @param message The message to print.
		 */
		void printMessage(const string & message);
	
	public:

		/**
		 * @brief Set the constraint policy for this optimizer.
		 *
		 * @param constraintPolicy The constraint policy.
		 */
		void setConstraintPolicy(const string & constraintPolicy);

		/**
		 * @brief Get the constraint policy for this optimizer.
		 *
		 * @return The constraint policy.
		 */
		string getConstraintPolicy() const;
		
		/**
		 * @brief Set the maximum number of function evaluation to perform during optimization.
		 *
		 * @param max The maximum number of evaluations to perform.
		 */
		void setMaximumNumberOfEvaluations(int max);
	
	public:
	
		static string CONSTRAINTS_AUTO;
		static string CONSTRAINTS_IGNORE;
		static string CONSTRAINTS_KEEP;
};

#endif	//_ABSTRACTOPTIMIZER_H_
