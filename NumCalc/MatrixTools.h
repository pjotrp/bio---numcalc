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
		if(ncA != nrB) throw DimensionException("MatrixTools::pow(). nrows B != ncols A.", nrB, ncA); 
		MatrixB C(nrA, ncB);
		mult(A, B, C);
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

};


#endif	//_MATRIXTOOLS_H_
