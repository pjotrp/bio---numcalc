//
// File: Matrix.h
// Created by: Julien Dutheil
// Created on: Tue Apr 07 11:58 2004
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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
	
		virtual const Scalar & operator()(unsigned int i, unsigned int j) const = 0;
		virtual       Scalar & operator()(unsigned int i, unsigned int j) = 0;
		virtual unsigned int nRows() const = 0;
		virtual unsigned int nCols() const = 0;
		virtual vector<Scalar> row(unsigned int i) const = 0;
		virtual vector<Scalar> col(unsigned int j) const = 0;
		virtual void resize(unsigned int nRows, unsigned int nCols) = 0;
};

template<class Scalar>
class RowMatrix : public Matrix<Scalar>, public vector< vector<Scalar> > {

	public:
		RowMatrix() {}

		RowMatrix(unsigned int nRow, unsigned int nCol): vector< vector<Scalar> >(nRow)
		{
			for(unsigned int i = 0; i < nRow; i++) {
				vector< vector<Scalar> >::operator[](i).resize(nCol);
			}
		}

		RowMatrix(const Matrix<Scalar> & m) 
		{
			resize(m.nRows(), m.nCols());
			for(unsigned int i = 0; i < m.nRows(); i++) {
				for(unsigned int j = 0; j < m.nCols(); j++) {
					operator()(i, j) = m(i, j);
				}
			}
		}

		~RowMatrix() {};

	public:
		
		Matrix<Scalar> & operator=(const Matrix<Scalar> & m)
		{
			vector< vector<Scalar> >::resize(m.nRows());	
			for(unsigned int i = 0; i < m.nRows(); i++) {
				vector< vector<Scalar> >::operator[](i).resize(m.nCols());
				for(unsigned int j = 0; j < m.nCols(); j++) {
					operator()(i,j) = m(i, j);
				}
			}
			return *this;
		}

		Clonable * clone() const { return new RowMatrix(* this); }
		
		const Scalar & operator()(unsigned int i, unsigned int j) const
		{
			return vector< vector<Scalar> >::operator[](i)[j];
		}
		
		Scalar & operator()(unsigned int i, unsigned int j)
		{
			return vector< vector<Scalar> >::operator[](i)[j];
		}
		
		unsigned int nRows() const { return vector< vector<Scalar> >::size(); }
		
		unsigned int nCols() const { return vector< vector<Scalar> >::size() == 0 ? 0 : vector< vector<Scalar> >::operator[](0).size(); }
		
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

		void resize(unsigned int nRows, unsigned int nCols)
		{
			vector< vector<double> >::resize(nRows);
			for(unsigned int i = 0; i < nRows; i++) {
				vector< vector<Scalar> >::operator[](i).resize(nCols);
			}
		}

};

#endif //_MATRIX_H_

