//
// File: MatrixTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Jan 19 16:42:25 2004
//

#ifndef _MATRIXTOOLS_H_
#define _MATRIXTOOLS_H_

#include "VectorTools.h"

//  From the MTL:
#include <mtl/mtl.h>
using namespace mtl;

namespace MatrixTools
{
	template<class Matrix>
	Matrix getId(unsigned int n)
	{
		Matrix id(n, n);
		for(unsigned int i = 0; i < n; i++) {
			for(unsigned int j = 0; j < n; j++) id(i, j) = (i == j) ? 1 : 0;
		}
		return id;
	}

	template<class Matrix, class T>
	Matrix diag(const vector<T> & d)
	{
		unsigned int n = d.size();
		Matrix diago(n, n);
		for(unsigned int i = 0; i < n; i++) {
			for(unsigned int j = 0; j < n; j++) diago(i, j) = (i == j) ? d[i] : 0;
		}
		return diago;
	}

	template<class Matrix, class T>
	vector<T> diag(const Matrix & M) throw (DimensionException)
	{
		unsigned int nc = M.ncols();
		unsigned int nr = M.nrows();
		if(nc != nr) throw DimensionException("MatrixTools::diag(). M must be a square matrix.", nr, nc); 
		vector<T> diago(nc);
		for(unsigned int i = 0; i < nc; i++) diago[i] = M(i, i);
		return diago;
	}

	template<class Matrix, class X>
	void fill(Matrix & M, X x)
	{
		for(unsigned int i = 0; i < M.nrows(); i++) {
			for(unsigned int j = 0; j < M.ncols(); j++) {
				M(i, j) = x;
			}
		}
	}

	template<class MatrixA, class MatrixB>
	MatrixB operator*(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		unsigned int ncA = A.ncols();
		unsigned int nrA = A.nrows();
		unsigned int nrB = B.nrows();
		unsigned int ncB = B.ncols();
		if(ncA != nrB) throw DimensionException("MatrixTools::operator*(). nrows B != ncols A.", nrB, ncA); 
		MatrixB C(nrA, ncB);
		mult(A, B, C);
		return C;
	}

	template<class MatrixA, class MatrixB>
	MatrixB operator+(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
	  unsigned int ncA = A.ncols();
		unsigned int nrA = A.nrows();
		unsigned int nrB = B.nrows();
		unsigned int ncB = B.ncols();
		if(ncA != ncB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of colums.", ncB, ncA); 
		if(nrA != nrB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of rows.", nrB, nrA); 
		MatrixB C(B);
		add(A, C);
		return C;
	}		
	
	template<class Matrix>
	Matrix operator-(const Matrix A)
	{
		Matrix B(A);
		scale(B, -1);
		return B;
	}

	template<class MatrixA, class MatrixB>
	MatrixB operator-(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
	  unsigned int ncA = A.ncols();
		unsigned int nrA = A.nrows();
		unsigned int nrB = B.nrows();
		unsigned int ncB = B.ncols();
		if(ncA != ncB) throw DimensionException("MatrixTools::operator-(). A and B must have the same number of colums.", ncB, ncA); 
		if(nrA != nrB) throw DimensionException("MatrixTools::operator-(). A and B must have the same number of rows.", nrB, nrA); 
		MatrixB C(-B);
		add(A, C);
		return C;
	}

	template<class Matrix>
	Matrix pow(const Matrix & m, int p) throw (DimensionException)
	{
		unsigned int n = m.nrows();
		if(n != m.ncols()) throw DimensionException("MatrixTools::pow(). nrows != ncols.", m.ncols(), m.nrows()); 
		if(p == 0) return getId<Matrix>(n);
		else {
			Matrix result(n, n);
			mult(pow<Matrix>(m, p - 1), m, result);
			return result;
		}
	}
	
	template<class Matrix, class Vector>
	Vector row(const Matrix & m, unsigned int rowIndex)
	{
		Vector r(m.ncols());
		for(unsigned int i = 0; i < m.ncols(); i++) r[i] = m(rowIndex, i);
		return(r);
	}

	template<class Matrix, class Vector>
	Vector col(const Matrix & m, unsigned int colIndex)
	{
		Vector c(m.nrows());
		for(unsigned int i = 0; i < m.nrows(); i++) c[i] = m(colIndex, i);
		return(c);
	}
	
};


#endif	//_MATRIXTOOLS_H_
