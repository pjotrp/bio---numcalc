//
// File: Matrix.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Apr 07 11:58 2004
//

/*
Copyright ou � ou Copr. Julien Dutheil, (17 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant � fournir des classes
pour le calcul num�rique.

Ce logiciel est r�gi par la licence CeCILL soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffus�e par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant 
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe � 
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement, 
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�. 

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accept� les
termes.
*/

/*
Copyright or � or Copr. Julien Dutheil, (November 17, 2004)

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


#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <vector>
using namespace std;

#include <Utils/Clonable.h>

/**
 * @brief The matrix template interface.
 */
template<class Scalar>
class Matrix: public Clonable {

	public:
		Matrix() {}
		virtual ~Matrix() {};

	public:
		template<class MatrixTemplate>
		Matrix<Scalar> & operator=(const MatrixTemplate & m) = 0;
		
		virtual const Scalar & operator()(unsigned int i, unsigned int j) const = 0;
		virtual       Scalar & operator()(unsigned int i, unsigned int j) = 0;
		virtual unsigned int nRows() const = 0;
		virtual unsigned int nCols() const = 0;
		virtual vector<Scalar> row(unsigned int i) const = 0;
		virtual vector<Scalar> col(unsigned int j) const = 0;
};

template<class Scalar>
class RowMatrix : public Matrix<Scalar>, public vector< vector<Scalar> > {

	public:
		RowMatrix() {}

		RowMatrix(unsigned int nRow, unsigned int nCol): vector< vector<Scalar> >(nRow)
		{
			for(unsigned int i = 0; i < nRow; i++) {
				operator[](i).resize(nCol);
			}
		}

		RowMatrix(const Matrix<Scalar> & m) 
		{
			resize(m.nRows());
			for(unsigned int i = 0; i < m.nRows(); i++) {
				operator[](i).resize(m.nCols());
				for(unsigned int j = 0; j < m.nCols(); j++) {
					operator()(i, j) = m(i, j);
				}
			}
		}

		~RowMatrix() {};

	public:
		
		template<class MatrixTemplate>
		Matrix<Scalar> & operator=(const MatrixTemplate & m)
		{
			resize(m.nRows());	
			for(unsigned int i = 0; i < m.nRows(); i++) {
				operator[](i).resize(m.nCols());
				for(unsigned int j = 0; j < m.nCols(); j++) {
					operator()(i,j) = m(i, j);
				}
			}
			return *this;
		}

		Clonable * clone() const { return new RowMatrix(* this); }
		
		const Scalar & operator()(unsigned int i, unsigned int j) const
		{
			return operator[](i)[j];
		}
		
		Scalar & operator()(unsigned int i, unsigned int j)
		{
			return operator[](i)[j];
		}
		
		unsigned int nRows() const { return size(); }
		
		unsigned int nCols() const { return size() == 0 ? 0 : operator[](0).size(); }
		
		vector<Scalar> row(unsigned int i) const
		{
			vector<Scalar> r(nCols());
			for(unsigned int j = 0; j < nCols(); j++) r[j] = operator()(i, j);
			return(r);
		}
		
		vector<Scalar> col(unsigned int j) const
		{
			vector<Scalar> c(nRows());
			for(unsigned int i = 0; i < nRows(); i++) c[i] = operator()(i, j);
			return(c);
		}

};

#endif //_MATRIX_H_
