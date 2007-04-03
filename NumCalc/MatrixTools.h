//
// File: MatrixTools.h
// Created by: Julien Dutheil
// Created on: Mon Jan 19 16:42:25 2004
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

/**
 * @brief Functions dealing with matrices.
 */
class MatrixTools
{
	public:
		MatrixTools() {}
		~MatrixTools() {}

	public:

		/**
		 * @return A copy of the given matrix.
		 * @param A original matrix.
		 */
		template<class Matrix>
		static Matrix copy(const Matrix & A)
		{
			return Matrix(A);
		}
	
		/**
		 * @return A identity matrix of size n.
		 * @param n the size of the matrix.
		 */
		template<class Matrix>
		static Matrix getId(unsigned int n)
		{
			Matrix id(n, n);
			for(unsigned int i = 0; i < n; i++)
      {
				for(unsigned int j = 0; j < n; j++) id(i, j) = (i == j) ? 1 : 0;
			}
			return id;
		}

		/**
		 * @return A diagonal matrix with diagonal elements taken from a vector.
		 * @param d A vector of diagonal elements.
		 */
		template<class Matrix, class T>
		static Matrix diag(const vector<T> & d)
		{
			unsigned int n = d.size();
			Matrix diago(n, n);
			for(unsigned int i = 0; i < n; i++)
      {
				for(unsigned int j = 0; j < n; j++) diago(i, j) = (i == j) ? d[i] : 0;
			}
			return diago;
		}

		/**
		 * @return The diagonal elements of a square matrix as a vector.
		 * @param M The matrix.
		 * @throw DimensionException If M is not a square matrix.
		 */
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

		/**
		 * @brief Set all elements in M to value x.
		 * @param M A matrix.
		 * @param x The value to use.
		 */
		template<class Matrix, class Scalar>
		static void fill(Matrix & M, Scalar x)
		{
			for(unsigned int i = 0; i < M.nRows(); i++)
      {
				for(unsigned int j = 0; j < M.nCols(); j++)
        {
					M(i, j) = x;
				}
			}
		}

		/**
		 * @brief Multiply all elements of a matrix by a given value, and add a constant.
		 *
		 * Performs \f$\forall i \forall j m_{i,j} = a.m_{i,j}+b\f$.
		 * 
		 * @param X A matrix.
		 * @param a Multiplicator.
		 * @param b Constant.
		 */
		template<class Matrix, class Scalar>
		static void scale(Matrix & X, Scalar a, Scalar b = 0)
    {
			for(unsigned int i = 0; i < X.nRows(); i++)
      {
				for(unsigned int j = 0; j < X.nCols(); j++)
        {
					X(i,j) = a * X(i, j) + b;
				}
			}
		}

		/**
		 * @return The dot product of two matrices.
		 * @param A First matrix.
		 * @param B Second matrix.
		 */
		template<class MatrixA, class MatrixB>
		static MatrixB mult(const MatrixA & A, const MatrixB & B) throw (DimensionException)
		{
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
			if(ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA); 
			MatrixB C(nrA, ncB);
			for(unsigned int i = 0; i < nrA; i++)
      {
				for(unsigned int j = 0; j < ncB; j++)
        {
					C(i, j) = 0;
					for(unsigned int k = 0; k < ncA; k++)
          {
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
		 * @param A The first matrix.
		 * @param D The diagonal matrix (only diagonal elements in a vector)
		 * @param B The second matrix.
		 * @return The result matrix.
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
			for(unsigned int i = 0; i < nrA; i++)
      {
				for(unsigned int j = 0; j < ncB; j++)
        {
					C(i, j) = 0;
					for(unsigned int k = 0; k < ncA; k++)
          {
						C(i, j) += A(i, k) * B(k, j) * D[k];
					}
				}
			}
			return C;
		}

		/**
		 * @brief Add matrix B to matrix A.
		 *
		 * @param A Matrix A
		 * @param B Matrix B
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
			for(unsigned int i = 0; i < A.nRows(); i++)
      {
				for(unsigned int j = 0; j < A.nCols(); j++)
        {
					A(i, j) += B(i, j);
				}
			}
		}		
	
		/**
		 * @return \f$\prod_{i=1}^p m\f$.
		 * @param m The matrix.
		 * @param p The number of multiplications.
		 * If p = 0, sends the identity matrix.
		 * @throw DimensionException If m is not a square matrix.
		 */
		template<class Matrix>
		static Matrix pow(const Matrix & m, unsigned int p) throw (DimensionException)
		{
			unsigned int n = m.nRows();
			if(n != m.nCols()) throw DimensionException("MatrixTools::pow(). nrows != ncols.", m.nCols(), m.nRows()); 
			if(p == 0) return getId<Matrix>(n);
			else
      {
				Matrix result(n, n);
				result = mult(m, pow<Matrix>(m, p - 1));
				return result;
			}
		}
	
		/**
		 * @return The position of the maximum value in the matrix.
		 * @param m The matrix.
		 */
		template<class Matrix>
		static vector<unsigned int> whichmax(const Matrix & m)
		{
			unsigned int nrows = m.nRows();
			unsigned int ncols = m.nCols();
			vector<unsigned int> pos(2);
			unsigned int imax = 0;
			unsigned int jmax = 0;
			double currentMax = log(0.);
			for(unsigned int i = 0; i < nrows; i++)
      {
				for(unsigned int j = 0; j < ncols; j++)
        {
					double currentValue = m(i, j);
					//cout << currentValue << "\t" << (currentValue > currentMax) << endl;
					if(currentValue > currentMax)
          {
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

		/**
		 * @return The position of the minimum value in the matrix.
		 * @param m The matrix.
		 */
		template<class Matrix>
		static vector<unsigned int> whichmin(const Matrix & m)
		{
			unsigned int nrows = m.nRows();
			unsigned int ncols = m.nCols();
			vector<unsigned int> pos(2);
			unsigned int imin = 0;
			unsigned int jmin = 0;
			double currentMin = -log(0.);
			for(unsigned int i = 0; i < nrows; i++)
      {
				for(unsigned int j = 0; j < ncols; j++)
        {
					double currentValue = m(i, j);
					if(currentValue < currentMin)
          {
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

		/**
		 * @brief Print a matrix to a stream.
		 * 
		 * @param m The matrix to print.
		 * @param out The stream to use.
		 */
		template<class Matrix>
		static void print(const Matrix & m, ostream & out = cout)
		{
			out << m.nRows() << "x" << m.nCols() << endl;
			out << "[" << endl;
			for(unsigned int i = 0; i < m.nRows(); i++)
      {
				out << "[";
				for(unsigned int j = 0; j < m.nCols() - 1; j++)
        {
					out << m(i, j) << ", ";
				}
				if(m.nCols() > 0) out << m(i, m.nCols() - 1) << "]" << endl;
			}
			out << "]" << endl;
		}
		
		/**
		 * @brief Print a vector to a stream.
		 * 
		 * @param v The vector to print.
		 * @param out The stream to use.
		 */
		template<class Real>
		static void print(const vector<Real> & v, ostream & out = cout)
		{
			out << v.size() << endl;
			out << "[";
			for(unsigned int i = 0; i < v.size() - 1; i++)
      {
				out << v[i] << ", ";
			}
			if(v.size() > 0) out << v[v.size() - 1];
			out << "]" << endl;
		}

		/**
		 * @return True if the matrix is a square matrix.
		 * @param A A matrix.
		 */
		template<class Matrix>
		static bool isSquare(const Matrix & A) { return A.nRows() == A.nCols(); }

		/**
		 * @return The inverse matrix of A.
		 * @param A The matrix to inverse.
		 * @throw DimensionException If A is not a square matrix.
		 */
		template<class Real>
		static RowMatrix<Real> inv(const Matrix<Real> & A) throw (DimensionException)
		{
			if(! isSquare(A)) throw DimensionException("MatrixTools::inv(). Matrix A is not a square matrix.", A.nRows(), A.nCols());
			LUDecomposition<Real> lu(A);
			RowMatrix<Real> I = getId<RowMatrix<Real> >(A.nRows());
			return lu.solve(I);
		}

		/**
		 * @return The transposition of A.
		 * @param A The matrix.
		 */
		template<class Real>
		static RowMatrix<Real> transpose(const Matrix<Real> & A)
		{
			RowMatrix<Real> M(A.nCols(), A.nRows());
			for(unsigned int i = 0; i < A.nCols(); i++)
      {
				for(unsigned int j = 0; j < A.nRows(); j++)
        {
					M(i, j) = A(j, i);
				}
			}
			return M;
		}
	
    /**
     * @brief Compute the Kronecker product of two row matrices.
     *
     * @param A The first row matrix.
     * @param B The second row matrix.
     * @return The product \f$A \otimes B\f$.
     */
    template<class Scalar>
    static RowMatrix<Scalar> kroneckerMult(const RowMatrix<Scalar> & A, const RowMatrix<Scalar> & B)
    {
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
      RowMatrix<Scalar> C(nrA*nrB, ncA*ncB);
      for(unsigned int ia = 0; ia < nrA; ia++)
      {
        for(unsigned int ja = 0; ja < ncA; ja++)
        {
          Scalar aij = A(ia, ja);
          for(unsigned int ib = 0; ib < nrB; ib++)
          {
            for(unsigned int jb = 0; jb < ncB; jb++)
            {
              C(ia*nrB+ib,ja*ncB+jb) = aij*B(ib,jb);
            }
          }
        }
      }
      return C;
    }

    /**
     * @brief Compute the Kronecker sum of two row matrices.
     *
     * @param A The first row matrix.
     * @param B The second row matrix.
     * @return The product \f$A \oplus B\f$.
     */
    template<class Scalar>
    static RowMatrix<Scalar> kroneckerSum(const RowMatrix<Scalar> & A, const RowMatrix<Scalar> & B)
    {
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
      RowMatrix<Scalar> C(nrA+nrB,ncA+ncB);
      for(unsigned int ia = 0; ia < nrA; ia++)
      {
        for(unsigned int ja = 0; ja < ncA; ja++)
        {
          C(ia,ja) = A(ia,ja);
        }
      }
      for(unsigned int ib = 0; ib < nrB; ib++)
      {
        for(unsigned int jb = 0; jb < nrB; jb++)
        {
          C(nrA+ib,ncA+jb) = B(ib,jb);
        }
      }
      return C;
    }

    /**
     * @brief Compute the Kronecker sum of n row matrices.
     *
     * @param A A vector of row matrices of any size.
     * @return The product \f$\bigoplus_i A_i\f$.
     */
    template<class Scalar>
    static RowMatrix<Scalar> kroneckerSum(const vector< RowMatrix<Scalar> > & A)
    {
			unsigned int nr = 0;
			unsigned int nc = 0;
      for(unsigned int k = 0; k < A.size(); k++)
      {
        nr += A[k].nRows();
        nc += A[k].nCols();
      }
      RowMatrix<Scalar> C(nr,nc);
      unsigned int rk = 0; //Row counter
      unsigned int ck = 0; //Col counter
      for(unsigned int k = 0; k < A.size(); k++)
      {
        const RowMatrix<Scalar> * Ak = &A[k];
        for(unsigned int i = 0; i < Ak->nRows(); i++)
        {
          for(unsigned int j = 0; j < Ak->nCols(); j++)
          {
            C(rk+i,ck+j) = (*Ak)(i,j);
          }
        }
        rk+=Ak->nRows();
        ck+=Ak->nCols();
      }
      return C;
    }

	
};

/* DEPRECATED 
namespace MatrixOperators {
	
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
*/

#endif	//_MATRIXTOOLS_H_
