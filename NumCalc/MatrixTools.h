//
// File: MatrixTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Jan 19 16:42:25 2004
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



#ifndef _MATRIXTOOLS_H_
#define _MATRIXTOOLS_H_

#include "VectorTools.h"
#include "Matrix.h"
#include "LUDecomposition.h"

#include <cstdio>
#include <iostream>
using namespace std;

class MatrixTools
{
	public:
		MatrixTools() {}
		~MatrixTools() {}

	public:

		template<class Matrix>
		static Matrix copy(const Matrix & A)
		{
			return Matrix(A);
		}
	
		template<class Matrix>
		static Matrix getId(unsigned int n)
		{
			Matrix id(n, n);
			for(unsigned int i = 0; i < n; i++) {
				for(unsigned int j = 0; j < n; j++) id(i, j) = (i == j) ? 1 : 0;
			}
			return id;
		}

		template<class Matrix, class T>
		static Matrix diag(const vector<T> & d)
		{
			unsigned int n = d.size();
			Matrix diago(n, n);
			for(unsigned int i = 0; i < n; i++) {
				for(unsigned int j = 0; j < n; j++) diago(i, j) = (i == j) ? d[i] : 0;
			}
			return diago;
		}

		template<class Matrix, class T>
		static vector<T> diag(const Matrix & M) throw (DimensionException)
		{
			unsigned int nc = M.nCols();
			unsigned int nr = M.nRows();
			if(nc != nr) throw DimensionException("MatrixTools::diag(). M must be a square matrix.", nr, nc); 
			vector<T> diago(nc);
			for(unsigned int i = 0; i < nc; i++) diago[i] = M(i, i);
			return diago;
		}

		template<class Matrix, class Scalar>
		static void fill(Matrix & M, Scalar x)
		{
			for(unsigned int i = 0; i < M.nRows(); i++) {
				for(unsigned int j = 0; j < M.nCols(); j++) {
					M(i, j) = x;
				}
			}
		}

		template<class Matrix, class Scalar>
		static void scale(Matrix & X, Scalar a, Scalar b = 0) {
			for(unsigned int i = 0; i < X.nRows(); i++) {
				for(unsigned int j = 0; j < X.nCols(); j++) {
					X(i,j) = a * X(i, j) + b;
				}
			}
		}

		template<class MatrixA, class MatrixB>
		static MatrixB mult(const MatrixA & A, const MatrixB & B) throw (DimensionException)
		{
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
			if(ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA); 
			MatrixB C(nrA, ncB);
			for(unsigned int i = 0; i < nrA; i++) {
				for(unsigned int j = 0; j < ncB; j++) {
					C(i, j) = 0;
					for(unsigned int k = 0; k < ncA; k++) {
						C(i, j) += A(i, k) * B(k, j);
					}
				}
			}
			return C;
		}

		/**
		 * @brief Compute A . D . B where D is a diagonal matrix in O(n^3).
		 *
		 * Since D is a diagonal matrix, this function is more efficient than doing
		 * mult(mult(A, diag(D)), B), which involves two 0(n^3) operations.
		 *
		 * @param MatrixA A The first matrix.
		 * @param vector<Scalar> D The diagonal matrix (only diagonal elements in a vector)
		 * @param MatrixB B The second matrix.
		 * @return MatrixB The result matrix.
		 * @throw DimensionException If matrices have not the appropriate size.
		 */
		template<class MatrixA, class Scalar, class MatrixB>
		static MatrixB mult(const MatrixA & A, const vector<Scalar> & D, const MatrixB & B) throw (DimensionException)
		{
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
			if(ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA); 
			if(ncA != D.size()) throw DimensionException("MatrixTools::mult(). Vector size is not eual to matrix size.", D.size(), ncA); 
			MatrixB C(nrA, ncB);
			for(unsigned int i = 0; i < nrA; i++) {
				for(unsigned int j = 0; j < ncB; j++) {
					C(i, j) = 0;
					for(unsigned int k = 0; k < ncA; k++) {
						C(i, j) += A(i, k) * B(k, j) * D[k];
					}
				}
			}
			return C;
		}

		/**
		 * @brief Add matrix B to matrix A.
		 *
		 * @param MatrixA A Matrix A
		 * @param MatrixB B Matrix B
		 * @throw DimensionException If A and B have note the same size.
		 */
		template<class MatrixA, class MatrixB>
		static void add(MatrixA & A, const MatrixB & B) throw (DimensionException)
		{
	 		unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
			if(ncA != ncB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of colums.", ncB, ncA); 
			if(nrA != nrB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of rows.", nrB, nrA); 
			for(unsigned int i = 0; i < A.nRows(); i++) {
				for(unsigned int j = 0; j < A.nCols(); j++) {
					A(i, j) += B(i, j);
				}
			}
		}		
	
		template<class Matrix>
		static Matrix pow(const Matrix & m, int p) throw (DimensionException)
		{
			unsigned int n = m.nRows();
			if(n != m.nCols()) throw DimensionException("MatrixTools::pow(). nrows != ncols.", m.nCols(), m.nRows()); 
			if(p == 0) return getId<Matrix>(n);
			else {
				Matrix result(n, n);
				result = m * pow<Matrix>(m, p - 1);
				return result;
			}
		}
	
		template<class Matrix>
		static vector<unsigned int> posmax(const Matrix & m)
		{
			unsigned int nrows = m.nRows();
			unsigned int ncols = m.nCols();
			vector<unsigned int> pos(2);
			unsigned int imax = 0;
			unsigned int jmax = 0;
			double currentMax = log(0.);
			for(unsigned int i = 0; i < nrows; i++) {
				for(unsigned int j = 0; j < ncols; j++) {
					double currentValue = m(i, j);
					//cout << currentValue << "\t" << (currentValue > currentMax) << endl;
					if(currentValue > currentMax) {
						imax = i;
						jmax = j;
						currentMax = currentValue;
					}
				}
			}
			pos[0] = imax;
			pos[1] = jmax;
			return pos;
		}

		template<class Matrix>
		static vector<unsigned int> posmin(const Matrix & m)
		{
			unsigned int nrows = m.nRows();
			unsigned int ncols = m.nCols();
			vector<unsigned int> pos(2);
			unsigned int imin = 0;
			unsigned int jmin = 0;
			double currentMin = -log(0.);
			for(unsigned int i = 0; i < nrows; i++) {
				for(unsigned int j = 0; j < ncols; j++) {
					double currentValue = m(i, j);
					if(currentValue < currentMin) {
						imin = i;
						jmin = j;
						currentMin = currentValue;
					}
				}
			}
			pos[0] = imin;
			pos[1] = jmin;
			return pos;
		}

		template<class Matrix>
		static void print(const Matrix & m, ostream & out = cout)
		{
			out << m.nRows() << "x" << m.nCols() << endl;
			out << "[" << endl;
			for(unsigned int i = 0; i < m.nRows(); i++) {
				out << "[";
				for(unsigned int j = 0; j < m.nCols() - 1; j++) {
					out << m(i, j) << ", ";
				}
				if(m.nCols() > 0) out << m(i, m.nCols() - 1) << "]" << endl;
			}
			out << "]" << endl;
		}
		
		template<class Real>
		static void print(const vector<Real> & v, ostream & out = cout)
		{
			out << v.size() << endl;
			out << "[";
			for(unsigned int i = 0; i < v.size() - 1; i++) {
				out << v[i] << ", ";
			}
			if(v.size() > 0) out << v[v.size() - 1];
			out << "]" << endl;
		}

		template<class Matrix>
		static bool isSquare(const Matrix & A) { return A.nRows() == A.nCols(); }

		template<class Real>
		static RowMatrix<Real> inv(const Matrix<Real> & A)
		{
			if(! isSquare(A)) throw DimensionException("MatrixTools::inv(). Matrix A is not a square matrix.", A.nRows(), A.nCols());
			LUDecomposition<Real> lu(A);
			RowMatrix<Real> I = getId<RowMatrix<Real> >(A.nRows());
			return lu.solve(I);
		}
	
	
};

namespace MatrixOperators {
	
	template<class MatrixA, class MatrixB>
	MatrixB operator*(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		return MatrixTools::mult<MatrixA, MatrixB>(A, B);
	}
	
	template<class MatrixA, class MatrixB>
	MatrixA operator+(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		MatrixA C = A;
		MatrixTools::add<MatrixA, MatrixB>(C, B);
		return C;
	}

	template<class MatrixA, class MatrixB>
	MatrixA operator+=(MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		MatrixTools::add<MatrixA, MatrixB>(A, B);
		return A;
	}

	template<class Matrix>
	Matrix operator-(const Matrix A)
	{
		Matrix B(A.nRows(), A.nCols());
		for(unsigned int i = 0; i < B.nRows(); i++) {
			for(unsigned int j = 0; j < B.nCols(); j++) {
				B(i, j) = -A(i, j);
			}
		}
		return B;
	}

	template<class MatrixA, class MatrixB>
	MatrixA operator-(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
//	  unsigned int ncA = A.nCols();
//		unsigned int nrA = A.nRows();
//		unsigned int nrB = B.nRows();
//		unsigned int ncB = B.nCols();
//		if(ncA != ncB) throw DimensionException("MatrixTools::operator-(). A and B must have the same number of colums.", ncB, ncA); 
//		if(nrA != nrB) throw DimensionException("MatrixTools::operator-(). A and B must have the same number of rows.", nrB, nrA); 
//		MatrixB C(A.nRows(), A.nCols());
//		for(unsigned int i = 0; i < A.nRows(); i++) {
//			for(unsigned int j = 0; j < A.nCols(); j++) {
//				C(i, j) = A(i, j) - B(i, j);
//			}
//		}
//		return C;
		MatrixA C = A;
		MatrixTools::add<MatrixA, MatrixB>(C, -B);
		return C;
	}
	
	template<class MatrixA, class MatrixB>
	MatrixA operator-=(MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		MatrixTools::add<MatrixA, MatrixB>(A, -B);
		return A;
	}

};

#endif	//_MATRIXTOOLS_H_
