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
MatrixTools.h\
Parameter.h\
Domain.h\
MultiDimensionsOptimizationTools.h\
ParameterList.h\
RandomFactory.h\
Uniform01K.h\
Uniform01QD.h\
Uniform01WH.h

CPP_FILES = \
AbstractOptimizer.cpp\
Functions.cpp\
Optimizer.cpp\
AutoParameter.cpp\
GammaDiscreteDistribution.cpp\
Parameter.cpp\
BrentOneDimension.cpp\
GoldenSectionSearch.cpp\
ParameterExceptions.cpp\
ConstantDistribution.cpp\
IntervalData.cpp\
ParameterList.cpp\
Constraints.cpp\
MultiDimensionsOptimizationTools.cpp\
Parametrizable.cpp\
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
VectorExceptions.cpp\
Uniform01K.cpp\
Uniform01QD.cpp\
Uniform01WH.cpp

lib_LTLIBRARIES = libnumcalc.la
libnumcalc_la_SOURCES = $(H_FILES) $(CPP_FILES)
libnumcalc_la_LDFLAGS = -version-info 1:0:0

library_includedir=$(includedir)/NumCalc
library_include_HEADERS = $(H_FILES)
