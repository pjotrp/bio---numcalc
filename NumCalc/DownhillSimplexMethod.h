//
// File: DownhillSimplexMethod.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov  4 17:10:05 2003
//

/*
Copyright ou © ou Copr. CNRS, (17 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

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

/**
 * @brief This implements the Downhill Simplex method in multidimensions.
 * The code is an adaptation of the one discribed in:
 * <pre>
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * (ISBN 0-521-43108-5)
 * </pre>
 */
class DownhillSimplexMethod: public AbstractOptimizer
{
	
	public:
		class DSMStopCondition: public AbstractOptimizationStopCondition
		{
			public:
				DSMStopCondition(DownhillSimplexMethod * dsm);
				~DSMStopCondition();
			
			public:
				void init() {}
				bool isToleranceReached() const;
		};
	
	friend class DSMStopCondition;
	
	protected:
		class Simplex: public vector<ParameterList> {
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
		int _iHighest, _iNextHighest, _iLowest;
	
	public: // constructor and destructor:

		/**
		 * @brief Build a new Downhill Simplex optimizer.
		 *
		 * @param function A pointer toward an object implementing the Optimizable interface.
		 */
		DownhillSimplexMethod(Function * function);
	
		virtual ~DownhillSimplexMethod();
	
	public: // The Optimizer interface:
		
		void init(const ParameterList & params) throw (Exception); //redefinition
	
		double step() throw (Exception);
	
		/**
		 * @brief Multidimensional minimization of the function _function by the
		 * downhill simplex method of Nelder and Mead.
		 *
		 * The simplex <i>p</i>[1..nDim+1][1..nDim]
		 * is input. Its <i>nDim+1</i> rows are nDim-dimensional vectors which are the vertices
		 * of the starting simplex. Also input is the vector <i>y</i>[1..nDim+1], whose components
		 * must be preinitialized to the values of _function evaluated at the <i>nDim + 1</i>
		 * vertices (rows). On output, <i>p</i> and <i>y</i> will have been reset to <i>nDim + 1</i>
		 * new points all within <i>fTol</i> of a minimum function value.
		 */
		double optimize() throw (Exception);
		double getFunctionValue() const throw (Exception);
	
	protected:
		
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
};


#endif	//_DOWNHILLSIMPLEXMETHOD_H_
