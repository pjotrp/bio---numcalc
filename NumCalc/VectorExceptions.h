#ifndef _VECTOREXCEPTIONS_H_
#define _VECTOREXCEPTIONS_H_

//From utils:
#include <Utils/Exceptions.h>

//From the STL/
#include <vector>
using namespace std;

/******************************************************************************/

template<class T> class VectorException : public Exception {

	protected:
		const vector<T> * _vect;
			
	public:
		// Class constructor
		VectorException(const char *   text, const vector<T> * vect = NULL) :
			Exception("VectorException: " + string(text)),
			_vect(vect) {};

		VectorException(const string & text, const vector<T> * vect = NULL) :
			Exception("VectorException: " + text),
			_vect(vect) {};

	
		// Class destructor
		~VectorException() throw () {};
		
	public:
		virtual const vector<T> * getVector() const { return _vect; }
};

/******************************************************************************/

template<class T> class EmptyVectorException : public VectorException<T> {

	public:
		// Class constructor
		EmptyVectorException(const char *   text, const vector<T> * vect = NULL) :
			VectorException<T>("EmptyVectorException: " + string(text), vect) {};
	
		EmptyVectorException(const string & text, const vector<T> * vect = NULL) :
			VectorException<T>("EmptyVectorException: " + text, vect) {};
	
		// Class destructor
		~EmptyVectorException() throw () {}
};

/******************************************************************************/

class DimensionException : public Exception {

	protected:
		int dimension;
		int correctDimension; 
			
	public:
		// Class constructor
		DimensionException(const char *   text, int dimension, int correctDimension);
		DimensionException(const string & text, int dimension, int correctDimension);
	
		// Class destructor
		~DimensionException() throw ();
	public:
		virtual int getDimension() const;
		virtual int getCorrectDimension() const;
};

/******************************************************************************/

#endif //_VECTOREXCEPTIONS_H_

