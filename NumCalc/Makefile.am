H_FILES = \
AbstractOptimizer.h\
DownhillSimplexMethod.h\
NumTools.h\
Parametrizable.h\
AutoParameter.h\
Functions.h\
OneDimensionOptimizationTools.h\
PowellMultiDimensions.h\
BrentOneDimension.h\
GammaDiscreteDistribution.h\
OptimizationStopCondition.h\
RandomTools.h\
ConstantDistribution.h\
GoldenSectionSearch.h\
Optimizer.h\
VectorTools.h\
VectorExceptions.h\
Constraints.h\
IntervalData.h\
ParameterExceptions.h\
DiscreteDistribution.h\
AbstractDiscreteDistribution.h\
SimpleDiscreteDistribution.h\
SimpleMultiDimensions.h\
Matrix.h\
MatrixTools.h\
EigenValue.h\
LUDecomposition.h\
Parameter.h\
Domain.h\
ParameterList.h\
RandomFactory.h\
Uniform01K.h\
Uniform01QD.h\
Uniform01WH.h\
DataTable.h\
DataTableExceptions.h\
parameters\
optimizers\
vectors\
matrix\
distributions\
random\
AbstractNumericalDerivative.h\
ThreePointsNumericalDerivative.h\
FivePointsNumericalDerivative.h

CPP_FILES = \
AbstractOptimizer.cpp\
AutoParameter.cpp\
GammaDiscreteDistribution.cpp\
Parameter.cpp\
BrentOneDimension.cpp\
GoldenSectionSearch.cpp\
ParameterExceptions.cpp\
ConstantDistribution.cpp\
IntervalData.cpp\
ParameterList.cpp\
AbstractDiscreteDistribution.cpp\
SimpleDiscreteDistribution.cpp\
SimpleMultiDimensions.cpp\
NumTools.cpp\
PowellMultiDimensions.cpp\
Domain.cpp\
OneDimensionOptimizationTools.cpp\
RandomTools.cpp\
DownhillSimplexMethod.cpp\
OptimizationStopCondition.cpp\
VectorTools.cpp\
Uniform01K.cpp\
Uniform01QD.cpp\
Uniform01WH.cpp\
DataTable.cpp\
ThreePointsNumericalDerivative.cpp\
FivePointsNumericalDerivative.cpp

lib_LTLIBRARIES = libbpp-numcalc.la
libbpp_numcalc_la_SOURCES = $(H_FILES) $(CPP_FILES)
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
libbpp_numcalc_la_LDFLAGS = -version-info 2:0:0

library_includedir=$(includedir)/NumCalc
library_include_HEADERS = $(H_FILES)

