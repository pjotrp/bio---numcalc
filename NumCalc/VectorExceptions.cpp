#include "VectorExceptions.h"

// From Utils:
#include <Utils/TextTools.h>
DimensionException::DimensionException(
	const char * text,
	int dimension,
	int correctDimension
) :
	Exception("DimensionException (found " +
		TextTools::toString(dimension) +
		", should be " +
		TextTools::toString(correctDimension) +
		string(text)),
	dimension(dimension),
	correctDimension(correctDimension) {};
		
DimensionException::DimensionException(
	const string & text,
	int dimension,
	int correctDimension
) :
	Exception("DimensionException (found " +
		TextTools::toString(dimension) +
		", should be " +
		TextTools::toString(correctDimension) +
		text),
	dimension(dimension),
	correctDimension(correctDimension) {};
		
DimensionException::~DimensionException() throw() {};

int DimensionException::getDimension() const { return dimension; }

int DimensionException::getCorrectDimension() const { return correctDimension; }

